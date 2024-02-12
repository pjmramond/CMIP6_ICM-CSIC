#!/bin/bash

# Input and output folder paths
input_folder="/Volumes/phalaropus/2023_SeaIcePenguinPopulations/ForcingData/01data/25km_withMembers/"
output_folder="/Volumes/phalaropus/2023_SeaIcePenguinPopulations/ForcingData/02regriddedData/"

crop_extent="-180,180,-90,-60"  # Adjust the extent as needed (lon_min, lon_max, lat_min, lat_max)
grid_file="grid_quarterDegree.txt"  # Specify the path to the grid file


# Function to process NetCDF files
process_nc_files() {
    local input_folder="$1"
    local output_folder="$2"

    # Recursively find all .nc files in the input folder while excluding hidden files
    find "$input_folder" -type f -name "*.nc" ! -name "._*" | while read -r input_file; do

        # Get the number of time steps for the current NetCDF file
        ntime=$(cdo ntime "$input_file")

        # Print information at the beginning of each iteration
        echo "Working on: $input_file, Time steps: $ntime"

        # Extract the relative path from the input folder
        relative_path="${input_file#$input_folder/}"

        # Construct the output folder path
        nc_output_folder="$output_folder/$relative_path"

        # Check if the directory already exists before attempting to create it
        if [ ! -d "$nc_output_folder" ]; then
            mkdir -p "$nc_output_folder"
        fi

        # Extract filename without extension
        filename=$(basename "$input_file" .nc)

        # Create a unique temporary regridded data file for each iteration
        temp_regridded_data="temp_regridded_data_${filename}.nc"

        # Regrid the input NetCDF file using remapcon with the specified grid_file
        cdo remapcon,"$grid_file" "$input_file" "$temp_regridded_data"

        # Crop and extract each time step as a separate NetCDF
        for ((i = 1; i <= ntime; i++)); do
            # Get the start date of the time step from the NetCDF file
            start_date=$(cdo showdate -seltimestep,$i "$temp_regridded_data" | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}')

            # Remove hyphens from the date
            start_date="${start_date//-/}"

            # Construct the final .nc file name with date at the end
            nc_file="$nc_output_folder/${filename}_$start_date.nc"

            # Crop the regridded NetCDF data to match the original extent and extract the time step
            cdo sellonlatbox,$crop_extent -seltimestep,$i "$temp_regridded_data" temp_cropped_data.nc

            # Write the cropped data to the hard drive
            cdo -copy temp_cropped_data.nc "$nc_file"

            echo "Processed and saved: $nc_file, Time step: $i of $ntime"
        done

        echo "Processed all time steps for $input_file"

        # Remove the unique temporary regridded data file
        rm "$temp_regridded_data"
    done

    # Remove temporary files
    rm temp_cropped_data.nc
}

# Call the function to process NetCDF files
process_nc_files "$input_folder" "$output_folder"
