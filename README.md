# Interpy My Bathy - Gap-Filling Interpolation Script

## Overview

This script is a graphical tool designed to fill gaps in bathymetric raster datasets (like GeoTIFF or BAG files). It intelligently identifies holes and missing data within your survey area and uses interpolation to create a continuous surface.

The primary workflow is a two-step process:
1.  **Generate a Preview Mask**: Based on your settings, the tool first shows you exactly which areas it plans to fill.
2.  **Run Interpolation**: If you approve the preview, the tool proceeds with the full, computationally intensive interpolation process.

This approach saves time by ensuring your parameters are correct before committing to the full run.

---

## How It Works

The script employs a sophisticated, multi-pass process to ensure efficient and high-quality gap filling, even on very large datasets.

### 1. Mask Generation

Before any interpolation happens, the script first defines exactly *where* to fill. It does this by creating a "final interpolation mask". You control this process with three key settings:

* **Bridging Distance**: This is the first step. The script expands the boundary of your existing data outwards by this distance. This is useful for connecting nearby survey lines or closing large inlets along the convex hull of your data.
* **Internal Gaps Only**: If this box is checked, the "Bridging" step is skipped. The script will only identify holes that are fully enclosed by existing data.
* **Max Gap Distance**: After the initial hull is created (with or without bridging), this parameter provides the final constraint. The script will only fill areas within the hull where a pixel is no more than this distance away from existing valid data. This prevents the tool from trying to invent data in the middle of extremely large, empty areas where interpolation would be unreliable.

### 2. Chunked Parallel Processing

To handle large raster files without running out of memory, the script doesn't try to process the entire file at once. Instead, it uses a chunking strategy.

* **Orientation Detection**: The script first analyzes the data to see if your survey lines are predominantly horizontal or vertical.
* **Smart Chunking**: It then slices the dataset into smaller, overlapping rectangular chunks perpendicular to the survey line orientation. This ensures that each chunk has the best possible data distribution for interpolation.
* **Parallel Processing**: These chunks are distributed across multiple CPU cores and are processed simultaneously (in parallel). This dramatically speeds up the interpolation process.
* **Blending**: The chunks overlap one another. The script smoothly blends the data in these overlapping zones to prevent any visible seams or edges in the final output.

### 3. Two-Pass Interpolation

The interpolation itself happens in two distinct passes to maximize quality:

* **Pass 1: Main Interpolation**: This is the primary run where the chunked parallel processing fills in all the gaps defined by the mask, using your chosen interpolation method (`linear`, `cubic`, or `nearest`).
* **Pass 2: Sliver Mop-Up**: Sometimes, the main interpolation can leave behind tiny, unfilled pixels or thin lines ("slivers"), especially at the seams between chunks. This second pass specifically detects these leftover artifacts and runs a localized, more focused interpolation to fill them in, ensuring a complete and clean final surface.

---

## How to Use the Tool

The graphical user interface (GUI) provides a straightforward way to control the interpolation process.

![image](bathynterpy.png)


### 1. File Selection
* **Input Raster File**: Click "Browse..." to select your input GeoTIFF or BAG file containing the bathymetry data with gaps.
* **Output Directory**: Click "Browse..." to choose the folder where the final interpolated GeoTIFF will be saved.

### 2. Parameter Configuration

These settings control the interpolation algorithm. The default values are a good starting point for moderately spaced survey lines.

* **Interpolation Method**:
    * `linear` (Default): A good balance of speed and smoothness. Connects known points with straight lines. Recommended for most cases.
    * `cubic`: Creates a smoother, more curved surface. It can be more computationally expensive and may sometimes produce values slightly outside the range of your source data.
    * `nearest`: The fastest method. It simply assigns the value of the nearest known pixel to a gap pixel. This results in a blocky, pixelated look and is best used for specific cases like categorical data, not elevation.

* **Max Gap Distance (m)**: The absolute maximum size of a gap that will be filled. A gap is defined by the distance from an empty pixel to the nearest valid data pixel.
    * *Default: 180.0*

* **Bridging Distance (m)**: The distance the tool will "reach out" from valid data to connect nearby areas before defining the outer boundary of the interpolation area. Set this to `0` if you only want to fill holes entirely enclosed by data.
    * *Default: 110.0*

* **Chunk Size (pixels)**: The size of the chunks the data is broken into for parallel processing. Larger values use more RAM but may be slightly more efficient with fewer chunks to manage. Smaller values use less RAM.
    * *Default: 512*

* **Overlap (pixels)**: How much the chunks should overlap to allow for smooth blending. This value should be large enough to cover any potential edge effects from the interpolation. It should generally not be changed unless you see seams in your output.
    * *Default: 200*

* **Fill Internal Gaps Only**: Check this box to disable the `Bridging Distance` functionality. The tool will only fill gaps that are fully surrounded by data and will not attempt to connect separate data areas.

### 3. Preview and Run

1.  Once your parameters are set, click the **"Preview Mask to Run Interpolation"** button.
2.  A plot window will appear, showing a preview of the mask. **White areas are the pixels that will be filled**. Black areas will be left untouched.
3.  Review the preview. If it looks correct, close the plot window and click **"Yes"** in the confirmation dialog box that appears.
4.  If the preview is not what you want, click **"No"**, adjust your `Max Gap Distance` or `Bridging Distance` parameters, and generate a new preview.
5.  After you confirm, the script will start the full interpolation process. You can monitor its progress in the console window. A message box will appear upon success or failure.

---

## Output Files

After a successful run, you will find the following file in your selected output directory:

* **`<input_filename>_<method>_interpolated.tif`**: This is your final, gap-filled GeoTIFF. If your input was a 2-band raster (e.g., elevation and uncertainty), this output will also have 2 bands. The uncertainty for newly interpolated areas will be set to a high value (99.0).
