# -*- coding: utf-8 -*-
"""
Gap-filling interpolation script for bathymetry datasets. Input is a geotiff or BAG; output is a geotiff.

Created on Tue Jun 17 15:02:00 2025
@author: Anthony.R.Klemm
"""

import numpy as np
import rasterio
import pyproj
from scipy.interpolate import griddata
from scipy.ndimage import binary_fill_holes, binary_dilation, binary_erosion, distance_transform_edt, label
import tkinter as tk
from tkinter import filedialog, Label, Entry, Button, StringVar, OptionMenu, messagebox, Checkbutton, BooleanVar
import os
import threading
import multiprocessing
from tqdm import tqdm
import matplotlib.pyplot as plt

def read_raster(file_path):
    print(f"[INFO] Reading raster from '{file_path}'...")
    with rasterio.open(file_path) as src:
        nodata_value = src.nodata if src.nodata is not None else 1000000.0
        data = src.read()
        profile = src.profile
        try:
            crs = pyproj.CRS(profile['crs'])
            if crs.is_compound:
                for sub_crs in crs.sub_crs_list:
                    if sub_crs.is_projected or sub_crs.is_geographic:
                        profile['crs'] = rasterio.crs.CRS.from_wkt(sub_crs.to_wkt())
                        break
        except Exception: pass
        profile.update(nodata=nodata_value)
        return data, nodata_value, profile

def create_interpolation_mask(data_band, profile, nodata_value,
                              max_gap_distance=50.0,
                              bridging_distance=15.0,
                              internal_only=False):
    print("[INFO] Creating final interpolation mask...")
    pixel_size = profile['transform'].a
    print(f"[DIAGNOSTIC] Pixel size: {pixel_size:.2f}m. Max gap: {max_gap_distance}m. Bridging: {bridging_distance}m.")

    original_footprint = (data_band != nodata_value) & (~np.isnan(data_band))

    if internal_only:
        print("[INFO] Internal Gaps Only mode selected. Skipping bridging.")
        hull_mask = binary_fill_holes(original_footprint)
    else:
        print(f"[INFO] Bridging mode selected. Bridging distance: {bridging_distance}m.")
        bridging_pixels = int(np.ceil(bridging_distance / pixel_size))
        structure = np.ones((3, 3), dtype=bool)
        bridged_mask = binary_dilation(original_footprint, structure=structure, iterations=bridging_pixels)
        filled_hull = binary_fill_holes(bridged_mask)
        hull_mask = binary_erosion(filled_hull, structure=structure, iterations=bridging_pixels)

    print("[INFO] Qualifying gaps to fill based on max distance...")
    distance_in_pixels = distance_transform_edt(~original_footprint)
    distance_mask = (distance_in_pixels * pixel_size) <= max_gap_distance

    calculated_fill_area = hull_mask & distance_mask
    print("[INFO] Re-combining with original data footprint to prevent edge trimming...")
    final_mask = calculated_fill_area | original_footprint

    pixels_to_fill = np.sum(final_mask) - np.sum(original_footprint)
    print(f"[DIAGNOSTIC] Preview mask generated. Found {pixels_to_fill:,} new pixels to fill.")
    return final_mask.astype(np.uint8)

def interpolate_chunk_worker(args):
    chunk_data, nodata_value, method, start_coord, end_coord, axis = args
    nodata_mask = (chunk_data == nodata_value) | np.isnan(chunk_data)
    if np.all(nodata_mask): return start_coord, end_coord, chunk_data
    valid_coords = np.array(np.nonzero(~nodata_mask)).T
    valid_values = chunk_data[~nodata_mask]
    rows, cols = np.indices(chunk_data.shape)
    interp_chunk = griddata(valid_coords, valid_values, (rows, cols), method=method, fill_value=nodata_value)
    return start_coord, end_coord, interp_chunk

def detect_orientation(data_mask):
    print("[INFO] Trying to detect survey line orientation...")
    horizontal_variance = np.var(data_mask, axis=1).sum()
    vertical_variance = np.var(data_mask, axis=0).sum()
    orientation = 'horizontal' if vertical_variance > horizontal_variance else 'vertical'
    print(f"[INFO] Detected predominantly {orientation} features.")
    return orientation

def chunked_interpolation(data_band, nodata_value, method='nearest', chunk_size=256, overlap=64):
    orientation = detect_orientation(data_band != nodata_value)
    if orientation == 'horizontal':
        return _chunked_interpolation_vertical(data_band, nodata_value, method, chunk_width=chunk_size, overlap=overlap)
    else:
        return _chunked_interpolation_horizontal(data_band, nodata_value, method, chunk_height=chunk_size, overlap=overlap)

def _run_parallel_interpolation(chunks_args_list):
    total_chunks = len(chunks_args_list)
    # num_processes = max(1, multiprocessing.cpu_count() - 1)
    num_processes = 8

    print(f"[INFO] Data split into {total_chunks} chunks. Starting a pool of {num_processes} worker processes.")
    results = []
    with multiprocessing.Pool(processes=num_processes) as pool:
        for result in tqdm(pool.imap(interpolate_chunk_worker, chunks_args_list), total=total_chunks, desc="Interpolating Chunks"):
            results.append(result)
    print("[INFO] Parallel processing complete.")
    return results

def _chunked_interpolation_horizontal(data_band, nodata_value, method, chunk_height, overlap):
    print(f"[INFO] Using HORIZONTAL chunking strategy (Chunk Height: {chunk_height}, Overlap: {overlap})...")
    height, width = data_band.shape
    chunks_args_list = []
    y = 0
    while y < height:
        start_y, end_y = y, min(y + chunk_height, height)
        chunks_args_list.append((data_band[start_y:end_y, :], nodata_value, method, start_y, end_y, 'y'))
        if end_y == height: break
        y += (chunk_height - overlap)
    interpolated_results = _run_parallel_interpolation(chunks_args_list)
    filled_data = np.full_like(data_band, nodata_value, dtype=np.float32)
    for i, (start_y, end_y, interp_chunk) in enumerate(sorted(interpolated_results)):
        write_start = overlap if i > 0 else 0
        filled_data[start_y + write_start:end_y, :] = interp_chunk[write_start:, :]
        if i > 0:
            pco = filled_data[start_y:start_y + overlap, :]
            cco = interp_chunk[:overlap, :]
            weights = np.linspace(1, 0, overlap, dtype=np.float32)[:, np.newaxis]
            blended = pco * weights + cco * (1 - weights)
            filled_data[start_y:start_y + overlap, :] = blended
    return filled_data

def _chunked_interpolation_vertical(data_band, nodata_value, method, chunk_width, overlap):
    print(f"[INFO] Using VERTICAL chunking strategy (Chunk Width: {chunk_width}, Overlap: {overlap})...")
    height, width = data_band.shape
    chunks_args_list = []
    x = 0
    while x < width:
        start_x, end_x = x, min(x + chunk_width, width)
        chunks_args_list.append((data_band[:, start_x:end_x], nodata_value, method, start_x, end_x, 'x'))
        if end_x == width: break
        x += (chunk_width - overlap)
    interpolated_results = _run_parallel_interpolation(chunks_args_list)
    filled_data = np.full_like(data_band, nodata_value, dtype=np.float32)
    for i, (start_x, end_x, interp_chunk) in enumerate(sorted(interpolated_results)):
        write_start = overlap if i > 0 else 0
        filled_data[:, start_x + write_start:end_x] = interp_chunk[:, write_start:]
        if i > 0:
            pco = filled_data[:, start_x:start_x + overlap]
            cco = interp_chunk[:, :overlap]
            weights = np.linspace(1, 0, overlap, dtype=np.float32)[np.newaxis, :]
            blended = pco * weights + cco * (1 - weights)
            filled_data[:, start_x:start_x + overlap] = blended
    return filled_data

def save_geotiff(file_path, data, profile, nodata_value):
    profile.update(driver='GTiff', dtype=rasterio.float32, count=data.shape[0], compress='lzw', nodata=nodata_value)
    with rasterio.open(file_path, 'w', **profile) as dst: dst.write(data.astype(np.float32))
    print(f"[INFO] GeoTIFF saved successfully to {file_path}")

def process_raster_threaded(input_file, output_dir, method, max_gap_distance, bridging_distance, chunk_size, overlap, internal_only, export_debug, interpolation_mask):
    try:
        print("\n[INFO] Starting full interpolation process with pre-computed mask...")
        data, nodata_value, profile = read_raster(input_file)
        num_bands = data.shape[0]
        base_name = os.path.basename(input_file)
        file_name, _ = os.path.splitext(base_name)
        output_file = os.path.join(output_dir, f"{file_name}_{method}_interpolated.tif")
        max_valid_elevation = 100.0
        nodata_tolerance = 0.001

        # --- Pass 1: Main Interpolation ---
        print("[INFO] Beginning Pass 1: Main parallel interpolation...")
        if num_bands == 1:
            band_to_interpolate = data[0]
        else: # num_bands == 2
            band_to_interpolate = data[0]
            uncertainty_band = data[1]

        filled_band = chunked_interpolation(band_to_interpolate, nodata_value, method, chunk_size, overlap)
        pass1_result = np.where(interpolation_mask == 1, filled_band, nodata_value)
        
        print("[INFO] Filtering high-value artifacts from interpolation result...")
        artifact_mask = (pass1_result > max_valid_elevation) & (np.abs(pass1_result - nodata_value) > nodata_tolerance)
        pass1_cleaned = np.where(artifact_mask, nodata_value, pass1_result)

        # --- Pass 2: Localized Re-interpolation for Slivers ---
        print("[INFO] Beginning Pass 2: Sliver detection and mop-up...")
        is_nodata_after_pass1 = (np.abs(pass1_cleaned - nodata_value) < nodata_tolerance) | np.isnan(pass1_cleaned)
        sliver_mask = (interpolation_mask == 1) & is_nodata_after_pass1

        if np.any(sliver_mask):
            labeled_slivers, num_blobs = label(sliver_mask)
            print(f"[DIAGNOSTIC] Found {np.sum(sliver_mask)} sliver pixels in {num_blobs} distinct blobs. Fixing them now...")
            
            final_band_to_process = pass1_cleaned.copy()
            
            for i in tqdm(range(1, num_blobs + 1), desc="Re-interpolating Slivers"):
                blob_mask = (labeled_slivers == i)
                rows, cols = np.where(blob_mask)
                margin = 20
                r_min, r_max = max(0, rows.min() - margin), min(final_band_to_process.shape[0], rows.max() + margin + 1)
                c_min, c_max = max(0, cols.min() - margin), min(final_band_to_process.shape[1], cols.max() + margin + 1)
                patch = final_band_to_process[r_min:r_max, c_min:c_max]
                patch_nodata_mask = (np.abs(patch - nodata_value) < nodata_tolerance) | np.isnan(patch)
                
                if not np.all(patch_nodata_mask):
                    patch_valid_coords = np.array(np.nonzero(~patch_nodata_mask)).T
                    patch_valid_values = patch[~patch_nodata_mask]
                    patch_rows, patch_cols = np.indices(patch.shape)

                    method_for_patch = method
                    # np.ptp (peak-to-peak) finds the range of values. If range is 0, all values are the same.
                    if len(patch_valid_coords) < 3 or np.ptp(patch_valid_coords[:, 0]) == 0 or np.ptp(patch_valid_coords[:, 1]) == 0:
                        if method != 'nearest':
                            print(f"[DIAGNOSTIC] Co-linear points detected for blob {i}. Temporarily falling back to 'nearest' for this patch.")
                        method_for_patch = 'nearest'

                    interp_patch = griddata(patch_valid_coords, patch_valid_values, (patch_rows, patch_cols), method=method_for_patch, fill_value=nodata_value)
                    local_blob_mask = blob_mask[r_min:r_max, c_min:c_max]
                    final_band_to_process[r_min:r_max, c_min:c_max][local_blob_mask] = interp_patch[local_blob_mask]

            final_elevation = final_band_to_process
        else:
            print("[INFO] No slivers found. Primary interpolation was complete.")
            final_elevation = pass1_cleaned
        
        # --- Final Assembly and Saving ---
        if num_bands == 1:
            save_geotiff(output_file, np.array([final_elevation]), profile, nodata_value)
        else: # num_bands == 2
            final_uncertainty = np.where(interpolation_mask == 1, uncertainty_band, 0.0)
            if np.any(sliver_mask): final_uncertainty[sliver_mask] = 99.0
            if np.any(artifact_mask): final_uncertainty = np.where(artifact_mask, 0.0, final_uncertainty)
            profile['count'] = 2
            save_geotiff(output_file, np.stack([final_elevation, final_uncertainty]), profile, nodata_value)

        # --- Optional Debug Export ---
        if export_debug:
            print("[INFO] Exporting debug rasters...")
            debug_profile = profile.copy()
            debug_profile.update(count=1, dtype=rasterio.uint8)
            preview_path = os.path.join(output_dir, f"{file_name}_preview_mask.tif")
            save_geotiff(preview_path, np.array([interpolation_mask]), debug_profile, 0)
            final_footprint_mask = (np.abs(final_elevation - nodata_value) > nodata_tolerance) & (~np.isnan(final_elevation))
            footprint_path = os.path.join(output_dir, f"{file_name}_final_footprint.tif")
            save_geotiff(footprint_path, np.array([final_footprint_mask.astype(np.uint8)]), debug_profile, 0)

        print("\n[SUCCESS] Process complete.")
        messagebox.showinfo("Success", f"Processing complete!\nOutput saved to:\n{output_file}")
    except Exception as e:
        print(f"[ERROR] An error occurred during the main processing: {e}")
        import traceback
        traceback.print_exc() # Print full traceback for debugging
        messagebox.showerror("Error", f"An error occurred during processing:\n{e}")

def select_input_file():
    file_path = filedialog.askopenfilename(title="Select Input Raster File", filetypes=[("Supported Rasters", "*.bag *.tif *.tiff"), ("All Files", "*.*")])
    if file_path: input_file_var.set(file_path)

def select_output_directory():
    directory = filedialog.askdirectory(title="Select Output Directory")
    if directory: output_dir_var.set(directory)

def run_preview_and_process():
    """Gets all parameters from GUI and starts the preview/process workflow."""
    input_file = input_file_var.get()
    output_dir = output_dir_var.get()
    if not input_file or not output_dir:
        messagebox.showwarning("Input Required", "Please select an input file and an output directory.")
        return
    try:
        max_gap = float(max_gap_dist_var.get())
        bridging = float(bridging_dist_var.get())
        chunk_size = int(chunk_size_var.get())
        overlap = int(overlap_var.get())
        method = method_var.get()
        internal_only = internal_only_var.get()
        export_debug = export_debug_var.get()
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter valid numbers for parameters.")
        return
    try:
        data, nodata_value, profile = read_raster(input_file)
        band_for_mask = data[0]
        band_profile = profile.copy()
        band_profile.update(count=1, height=band_for_mask.shape[0], width=band_for_mask.shape[1])
        interpolation_mask = create_interpolation_mask(band_for_mask, band_profile, nodata_value, max_gap, bridging, internal_only)
        plt.figure(figsize=(10, 8))
        plt.imshow(interpolation_mask, cmap='gray')
        plt.title("Interpolation Mask Preview\n(White areas will be filled)")
        plt.show(block=False)
    except Exception as e:
        messagebox.showerror("Preview Error", f"Could not generate preview:\n{e}")
        return
    proceed = messagebox.askyesno("Confirm Parameters", "The preview mask has been generated.\n\nDo you want to continue with the full interpolation?")
    if proceed:
        plt.close()
        threading.Thread(target=process_raster_threaded, args=(
            input_file, output_dir, method, max_gap, bridging, chunk_size, overlap, internal_only, export_debug, interpolation_mask
        ), daemon=True).start()
    else:
        plt.close()
        print("[INFO] User cancelled full processing.")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    root = tk.Tk()
    root.title("Interpy my bathy")
    
    input_file_var, output_dir_var = StringVar(), StringVar()
    method_var = StringVar(value='linear')
    max_gap_dist_var = StringVar(value='180.0')
    bridging_dist_var = StringVar(value='110.0')
    chunk_size_var = StringVar(value='512')
    overlap_var = StringVar(value='200')
    internal_only_var = BooleanVar(value=False)
    export_debug_var = BooleanVar(value=False)
    
    main_frame = tk.Frame(root, padx=10, pady=10)
    main_frame.pack(fill="both", expand=True)
    main_frame.columnconfigure(1, weight=1)

    row_index = 0
    Label(main_frame, text="Input Raster File:").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=input_file_var, width=60).grid(row=row_index, column=1, sticky="ew", pady=2)
    Button(main_frame, text="Browse...", command=select_input_file).grid(row=row_index, column=2, padx=5, pady=2)
    row_index += 1
    
    Label(main_frame, text="Output Directory:").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=output_dir_var, width=60).grid(row=row_index, column=1, sticky="ew", pady=2)
    Button(main_frame, text="Browse...", command=select_output_directory).grid(row=row_index, column=2, padx=5, pady=2)
    row_index += 1

    Label(main_frame, text="Interpolation Method:").grid(row=row_index, column=0, sticky="w", pady=2)
    OptionMenu(main_frame, method_var, 'nearest', 'linear', 'cubic').grid(row=row_index, column=1, sticky="w", pady=2)
    row_index += 1

    Label(main_frame, text="Max Gap Distance (m):").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=max_gap_dist_var).grid(row=row_index, column=1, sticky="w", pady=2)
    row_index += 1

    Label(main_frame, text="Bridging Distance (m):").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=bridging_dist_var).grid(row=row_index, column=1, sticky="w", pady=2)
    row_index += 1
    
    Label(main_frame, text="Chunk Size (pixels):").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=chunk_size_var).grid(row=row_index, column=1, sticky="w", pady=2)
    row_index += 1

    Label(main_frame, text="Overlap (pixels):").grid(row=row_index, column=0, sticky="w", pady=2)
    Entry(main_frame, textvariable=overlap_var).grid(row=row_index, column=1, sticky="w", pady=2)
    row_index += 1
    
    Checkbutton(main_frame, text="Fill Internal Gaps Only (No Bridging)", variable=internal_only_var).grid(row=row_index, column=1, sticky="w", padx=5, pady=2)
    row_index += 1
    
    #only needed for debugging
    # Checkbutton(main_frame, text="Export Debug Rasters", variable=export_debug_var).grid(row=row_index, column=1, sticky="w", padx=5, pady=2)
    # row_index += 1

    Button(main_frame, text="Preview Mask to Run Interpolation", command=run_preview_and_process, bg="lightblue").grid(row=row_index, column=1, columnspan=1, pady=20, sticky="w")
    
    root.mainloop()
