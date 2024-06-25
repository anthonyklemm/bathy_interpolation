import numpy as np
import rasterio
from scipy.ndimage import distance_transform_edt
from rasterio.warp import reproject, Resampling
import geopandas as gpd
from rasterio.features import shapes, geometry_mask
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
from skimage.morphology import binary_dilation, binary_erosion
from scipy.interpolate import griddata
from rasterio.transform import Affine
import threading
import tkinter as tk
from tkinter import filedialog, Label, Entry, Button, StringVar, OptionMenu
import multiprocessing

def read_raster(file_path):
    with rasterio.open(file_path) as src:
        nodata_value = src.nodata if src.nodata is not None else 1000000
        data = src.read()  # Read all bands
        profile = src.profile
        profile.update(nodata=nodata_value)
    return data, nodata_value, profile

def resample_raster(data, profile, scale_factor=0.1):
    """ Resample raster to a coarser resolution to simplify geometry creation. """
    new_width = int(data.shape[2] * scale_factor)
    new_height = int(data.shape[1] * scale_factor)
    resampled_data = np.empty((data.shape[0], new_height, new_width), dtype=np.float32)

    resampled_transform = profile['transform'] * Affine.scale(
        (data.shape[2] / new_width),
        (data.shape[1] / new_height)
    )

    for i in range(data.shape[0]):
        reproject(
            source=data[i],
            destination=resampled_data[i],
            src_transform=profile['transform'],
            src_crs=profile['crs'],
            dst_transform=resampled_transform,
            dst_crs=profile['crs'],
            resampling=Resampling.bilinear
        )

    return resampled_data, resampled_transform

def nearest_neighbor_interpolation(data, nodata_value):
    nodata_mask = data == nodata_value
    coords = np.indices(data.shape).reshape(2, -1).T
    valid_coords = np.array(np.nonzero(~nodata_mask)).T
    valid_values = data[~nodata_mask]

    filled_data = griddata(valid_coords, valid_values, coords, method='nearest', fill_value=nodata_value).reshape(data.shape)
    return filled_data

def bilinear_interpolation(data, nodata_value):
    nodata_mask = data == nodata_value
    coords = np.indices(data.shape).reshape(2, -1).T
    valid_coords = np.array(np.nonzero(~nodata_mask)).T
    valid_values = data[~nodata_mask]

    filled_data = griddata(valid_coords, valid_values, coords, method='linear', fill_value=nodata_value).reshape(data.shape)
    return filled_data

def cubic_interpolation(data, nodata_value):
    nodata_mask = data == nodata_value
    coords = np.indices(data.shape).reshape(2, -1).T
    valid_coords = np.array(np.nonzero(~nodata_mask)).T
    valid_values = data[~nodata_mask]

    filled_data = griddata(valid_coords, valid_values, coords, method='cubic', fill_value=nodata_value).reshape(data.shape)
    return filled_data

def apply_dilation_erosion(binary_array, dilation_iterations, erosion_iterations):
    for _ in range(dilation_iterations):
        binary_array = binary_dilation(binary_array)
    for _ in range(erosion_iterations):
        binary_array = binary_erosion(binary_array)
    return binary_array

def create_binary_mask_from_polygon(data, transform, nodata_value):
    binary_mask = (data != nodata_value).astype(np.uint8)
    binary_mask = apply_dilation_erosion(binary_mask, dilation_iterations=8, erosion_iterations=8)
    binary_mask = binary_mask.astype(np.uint8)
    polygons = [shape(geom) for geom, val in shapes(binary_mask, mask=binary_mask, transform=transform) if val]
    polygon = unary_union(polygons)
    return polygon, geometry_mask([mapping(polygon)], transform=transform, out_shape=(data.shape[1], data.shape[2]), invert=True)

def save_geotiff(file_path, data, profile):
    profile.update(dtype=rasterio.float32, count=data.shape[0], compress='lzw')
    with rasterio.open(file_path, 'w', **profile) as dst:
        for i in range(data.shape[0]):
            dst.write(data[i].astype(rasterio.float32), i + 1)

def chunk_interpolation(data, nodata_value, method='nearest', chunk_height=1500, overlap=600):
    if data.ndim == 2:  # Single band case
        data = np.expand_dims(data, axis=0)

    bands, height, width = data.shape
    filled_data = np.copy(data)

    # Calculate the correct number of chunks considering the overlap
    num_main_chunks_y = (height + chunk_height - 1) // chunk_height
    total_chunks = num_main_chunks_y + (num_main_chunks_y - 1)  # Main chunks + overlapping chunks

    chunk_count = 0

    # Process main chunks
    for y in range(0, height, chunk_height):
        chunk_count += 1
        print(f"Processing chunk {chunk_count} out of {total_chunks}")
        y_end = min(y + chunk_height, height)

        for band in range(bands):
            chunk = data[band, y:y_end, :]

            if np.all(chunk == nodata_value):
                # Skip this chunk as it contains only nodata values
                continue

            if method == 'nearest':
                filled_chunk = nearest_neighbor_interpolation(chunk, nodata_value)
            elif method == 'bilinear':
                filled_chunk = bilinear_interpolation(chunk, nodata_value)
            elif method == 'cubic':
                filled_chunk = cubic_interpolation(chunk, nodata_value)
            else:
                raise ValueError(f"Unsupported interpolation method: {method}")

            filled_data[band, y:y_end, :] = filled_chunk

    # Process overlapping chunks
    for y in range(chunk_height - overlap, height, chunk_height):
        chunk_count += 1
        print(f"Processing chunk {chunk_count} out of {total_chunks}")
        y_start = max(y - overlap, 0)
        y_end = min(y + overlap, height)

        for band in range(bands):
            chunk = data[band, y_start:y_end, :]

            if np.all(chunk == nodata_value):
                # Skip this chunk as it contains only nodata values
                continue

            if method == 'nearest':
                filled_chunk = nearest_neighbor_interpolation(chunk, nodata_value)
            elif method == 'bilinear':
                filled_chunk = bilinear_interpolation(chunk, nodata_value)
            elif method == 'cubic':
                filled_chunk = cubic_interpolation(chunk, nodata_value)
            else:
                raise ValueError(f"Unsupported interpolation method: {method}")

            # Determine the interior portion to update (excluding the edges)
            y_overlap_start = max(y, y + overlap // 2)
            y_overlap_end = min(y + overlap, y_end - overlap // 2)

            # Only update the interior portion in the filled_data array
            filled_data[band, y_overlap_start:y_overlap_end, :] = filled_chunk[
                y_overlap_start - y_start: y_overlap_end - y_start, :]

    return filled_data.squeeze()  # Remove the added dimension for single-band data


def process_band(input_file, output_dir, method, band_idx, nodata_value, profile):
    print(f"Processing band {band_idx + 1}...")
    data, _, _ = read_raster(input_file)
    data = data[band_idx]  # Extract the specific band
    data = np.expand_dims(data, axis=0)  # Add the band dimension back for consistency
    data_resampled, transform_resampled = resample_raster(data, profile)
    polygon, binary_mask_resampled = create_binary_mask_from_polygon(data_resampled, transform_resampled, nodata_value)
    binary_mask = np.empty((data.shape[1], data.shape[2]), dtype=np.uint8)
    reproject(
        source=binary_mask_resampled.astype(np.uint8),
        destination=binary_mask,
        src_transform=transform_resampled,
        src_crs=profile['crs'],
        dst_transform=profile['transform'],
        dst_crs=profile['crs'],
        resampling=Resampling.nearest
    )
    filled_data = chunk_interpolation(data, nodata_value, method=method)
    combined_data = np.where(binary_mask == 1, filled_data, nodata_value)
    return combined_data

def main(input_file, output_dir, method='nearest'):
    data, nodata_value, profile = read_raster(input_file)
    num_bands = data.shape[0]
    processed_bands = []

    for band_idx in range(num_bands):
        processed_band = process_band(input_file, output_dir, method, band_idx, nodata_value, profile)
        processed_bands.append(processed_band)

    combined_data = np.stack(processed_bands)
    output_file = f"{output_dir}/{input_file.split('/')[-1].split('.')[0]}_{method}.tif"
    print("Saving output GeoTIFF...")
    save_geotiff(output_file, combined_data, profile)
    print(f"Gap filling and clipping complete. The interpolated GeoTIFF has been saved as '{output_file}'.")

def select_input_file():
    file_path = filedialog.askopenfilename(
        title="Select Input File",
        filetypes=[("Raster Files", "*.bag *.tiff *.tif"), ("All Files", "*.*")]
    )
    if file_path:
        input_file_var.set(file_path)

def select_output_directory():
    directory = filedialog.askdirectory(title="Select Output Directory")
    if directory:
        output_dir_var.set(directory)

def run_interpolation():
    input_file = input_file_var.get()
    output_dir = output_dir_var.get()
    method = method_var.get()
    if input_file and output_dir:
        threading.Thread(target=main, args=(input_file, output_dir, method)).start()

root = tk.Tk()
root.title("Raster Interpolation GUI")

input_file_var = StringVar()
output_dir_var = StringVar()
method_var = StringVar(value='nearest')

Label(root, text="Input File:").grid(row=0, column=0, padx=10, pady=5)
Entry(root, textvariable=input_file_var, width=50).grid(row=0, column=1, padx=10, pady=5)
Button(root, text="Browse", command=select_input_file).grid(row=0, column=2, padx=10, pady=5)

Label(root, text="Output Directory:").grid(row=1, column=0, padx=10, pady=5)
Entry(root, textvariable=output_dir_var, width=50).grid(row=1, column=1, padx=10, pady=5)
Button(root, text="Browse", command=select_output_directory).grid(row=1, column=2, padx=10, pady=5)

Label(root, text="Interpolation Method:").grid(row=2, column=0, padx=10, pady=5)
OptionMenu(root, method_var, 'nearest', 'bilinear', 'cubic').grid(row=2, column=1, padx=10, pady=5)

Button(root, text="Run", command=run_interpolation).grid(row=3, column=1, pady=20)

root.mainloop()
