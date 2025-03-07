#!/bin/bash

# ======== CONFIGURATION SETTINGS ========
# Script version
VERSION="1.0.0"

# Flag to attempt recovery of technical sequences (I1/I2) when missing
# Set to "true" if you want to attempt to recover index files when they appear to be missing
RECOVER_TECHNICAL_READS=${RECOVER_TECHNICAL_READS:-false}

# Flag to attempt fixing malformed FASTQ files
# Set to "true" to attempt repair of problematic FASTQ files (default on)
FIX_FASTQ_FORMAT=${FIX_FASTQ_FORMAT:-true}

# Ensure script exits on error
set -e

# User-configurable settings (can be overridden by environment variables)
THREADS_PER_JOB=${THREADS_PER_JOB:-12}
PARALLEL_JOBS=${PARALLEL_JOBS:-4}
BASE_DIR="$(pwd)"
STATUS_DIR="${BASE_DIR}/.status"
# User should set this to their cellranger reference path
CELLRANGER_REF=${CELLRANGER_REF:-"/path/to/refdata-gex-GRCh38-2024-A"}
OUTPUT_DIR=${CUSTOM_OUTPUT_DIR:-"${BASE_DIR}"}
KEEP_FASTQS=${KEEP_FASTQS:-false}
DOWNLOAD_TIMEOUT=${DOWNLOAD_TIMEOUT:-3600}  # 1 hour timeout for downloads

# ======== UTILITY FUNCTIONS ========
# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check dependencies
check_dependencies() {
    log_message "Checking for required tools..."
    
    # Check bash version (need 4.0+ for associative arrays)
    if [ "${BASH_VERSINFO[0]}" -lt 4 ]; then
        log_message "Error: This script requires Bash 4.0 or higher"
        exit 1
    fi
    
    # Check for kingfisher
    if ! command -v kingfisher &> /dev/null; then
        log_message "Error: kingfisher is not installed or not in your PATH"
        log_message "Please install kingfisher with: pip install kingfisher-download"
        exit 1
    fi
    
    # Check for cellranger
    if ! command -v cellranger &> /dev/null; then
        log_message "Error: cellranger is not installed or not in your PATH"
        exit 1
    fi

    # Check for GNU awk (provides more features for FASTQ fixing)
    if ! command -v awk &> /dev/null; then
        log_message "Error: awk not found. Please install GNU awk."
        exit 1
    fi

    # Check for other essential tools
    for cmd in find grep sort wc curl wget gzip zcat sed stat mktemp; do
        if ! command -v $cmd &> /dev/null; then
            log_message "Error: $cmd not found. Please install it."
            exit 1
        fi
    done
    
    log_message "All required tools are available"
}

# Check input files
check_inputs() {
    if [ ! -f "SraRunTable.txt" ]; then
        log_message "Error: SraRunTable.txt not found in current directory"
        log_message "This file should be a tab-delimited file containing run metadata"
        log_message "You can download it from the NCBI SRA Run Selector"
        exit 1
    fi
    
    if [ ! -d "$CELLRANGER_REF" ]; then
        log_message "Error: Cellranger reference path does not exist: $CELLRANGER_REF"
        log_message "Please set CELLRANGER_REF to a valid directory"
        exit 1
    fi
}

# Function to validate FASTQ format
validate_fastq() {
    local fastq_file="$1"
    local check_lines=${2:-1000}  # Number of entries to check (default: 1000)
    
    log_message "Validating FASTQ format for $fastq_file..."
    
    # Check if the file exists and is not empty
    if [ ! -s "$fastq_file" ]; then
        log_message "Error: FASTQ file $fastq_file does not exist or is empty"
        return 1
    fi
    
    # Create a temporary file to store uncompressed data for checking
    local temp_file=$(mktemp)
    
    # Extract a sample of the file for validation
    zcat "$fastq_file" | head -n $((check_lines * 4)) > "$temp_file"
    
    # Check the basic FASTQ format
    local line_count=$(wc -l < "$temp_file")
    if [ $((line_count % 4)) -ne 0 ]; then
        log_message "Error: FASTQ file $fastq_file has $line_count lines, which is not divisible by 4"
        rm "$temp_file"
        return 1
    fi
    
    # Check the format of each entry
    local valid=true
    local line_num=1
    while [ "$line_num" -le "$line_count" ] && [ "$valid" = true ]; do
        # Line 1 should start with @
        local line1=$(sed -n "${line_num}p" "$temp_file")
        if [[ ! "$line1" =~ ^@ ]]; then
            log_message "Error: Line $line_num should start with @ in $fastq_file"
            valid=false
        fi
        
        # Line 3 should be a + sign (possibly followed by the same identifier)
        local line3=$(sed -n "$((line_num + 2))p" "$temp_file")
        if [[ ! "$line3" =~ ^\+ ]]; then
            log_message "Error: Line $((line_num + 2)) should start with + in $fastq_file"
            valid=false
        fi
        
        # Move to the next entry
        line_num=$((line_num + 4))
    done
    
    # Clean up
    rm "$temp_file"
    
    if [ "$valid" = true ]; then
        log_message "FASTQ validation passed for $fastq_file"
        return 0
    else
        log_message "FASTQ validation failed for $fastq_file"
        return 1
    fi
}

# Function to fix malformed FASTQ files
fix_fastq_format() {
    local fastq_file="$1"
    local fixed_file="${fastq_file}.fixed"
    
    log_message "Attempting to fix malformed FASTQ file: $fastq_file"
    
    # Create a temporary directory
    local temp_dir=$(mktemp -d)
    local temp_extracted="${temp_dir}/extracted.fastq"
    
    # Extract the file
    zcat "$fastq_file" > "$temp_extracted"
    
    # Check the file size (avoid processing huge files in memory)
    local file_size=$(stat -c %s "$temp_extracted")
    if [ "$file_size" -gt 10000000000 ]; then  # 10GB limit
        log_message "Error: File too large for in-memory fixing (${file_size} bytes)"
        rm -rf "$temp_dir"
        return 1
    fi
    
    # Count the number of @ header lines
    local header_count=$(grep -c "^@" "$temp_extracted")
    log_message "Found $header_count sequence headers"
    
    # Count the number of + separator lines
    local plus_count=$(grep -c "^+" "$temp_extracted")
    log_message "Found $plus_count separator lines"
    
    if [ "$header_count" -ne "$plus_count" ]; then
        log_message "Error: Mismatch between header count ($header_count) and separator count ($plus_count)"
        rm -rf "$temp_dir"
        return 1
    fi
    
    # Fix the file by ensuring every 4 lines follow FASTQ format
    awk '
    BEGIN { line_count = 0; }
    {
        line_count++;
        if (line_count % 4 == 1 && !($0 ~ /^@/)) {
            print "@" $0;  # Add @ if missing from header
        } else if (line_count % 4 == 3 && !($0 ~ /^\+/)) {
            print "+";     # Replace with simple + separator
        } else {
            print $0;      # Print unchanged
        }
    }' "$temp_extracted" > "${temp_dir}/fixed.fastq"
    
    # Compress the fixed file
    gzip -c "${temp_dir}/fixed.fastq" > "$fixed_file"
    
    # Validate the fixed file
    if validate_fastq "$fixed_file" 5000; then
        log_message "Successfully fixed and validated FASTQ file"
        mv "$fixed_file" "$fastq_file"  # Replace original with fixed
        rm -rf "$temp_dir"
        return 0
    else
        log_message "Failed to fix FASTQ file format"
        rm -rf "$temp_dir"
        rm -f "$fixed_file"
        return 1
    fi
}

# ======== SRA PROCESSING FUNCTIONS ========
# Function to check if FASTQ files already exist
check_fastq_files_exist() {
    local srr_accession="$1"
    local library_name="$2"
    
    # Check in the output directory for this library
    if [ -d "${OUTPUT_DIR}/${library_name}" ]; then
        local fastq_count=$(find "${OUTPUT_DIR}/${library_name}" -name "*${srr_accession}*.fastq.gz" 2>/dev/null | wc -l)
        if [ "$fastq_count" -gt 0 ]; then
            log_message "FASTQ files for $srr_accession already exist in ${OUTPUT_DIR}/${library_name}"
            return 0
        fi
    fi
    
    # Check in the SRR directory
    if [ -d "${srr_accession}" ]; then
        local fastq_count=$(find "${srr_accession}" -name "*.fastq.gz" 2>/dev/null | wc -l)
        if [ "$fastq_count" -gt 0 ]; then
            log_message "FASTQ files for $srr_accession already exist in ${srr_accession}/"
            return 0
        fi
    fi
    
    return 1
}

# Function to process a single SRA accession with kingfisher
process_sra_dir() {
    local srr_accession="$1"
    local library_name="$2"
    
    # Skip if already processed
    if [ -f "${STATUS_DIR}/${srr_accession}.processed" ]; then
        log_message "Skipping $srr_accession - already processed"
        return 0
    fi
    
    # Check if FASTQ files already exist
    if check_fastq_files_exist "$srr_accession" "$library_name"; then
        log_message "FASTQ files already exist for $srr_accession, marking as processed"
        touch "${STATUS_DIR}/${srr_accession}.processed"
        return 0
    fi
    
    log_message "Processing accession: $srr_accession for library: $library_name"
    
    # Create directory
    mkdir -p "$srr_accession"
    
    # Download with kingfisher
    log_message "Downloading $srr_accession with kingfisher..."
    
    # First attempt: Standard kingfisher download
    if ! timeout $DOWNLOAD_TIMEOUT kingfisher get -r "$srr_accession" \
                      -m ena-ftp aws-http prefetch \
                      --download-threads "$THREADS_PER_JOB" \
                      -t "$THREADS_PER_JOB" \
                      -f fastq.gz \
                      --output-directory "$srr_accession"; then
        log_message "Error: kingfisher download failed for $srr_accession"
        return 1
    fi
    
    # Check if download succeeded
    if [ "$(find "$srr_accession" -name "*.fastq.gz" | wc -l)" -eq 0 ]; then
        # Check if files are in a subdirectory
        if [ -d "$srr_accession/fastq" ] && [ "$(find "$srr_accession/fastq" -name "*.fastq.gz" | wc -l)" -gt 0 ]; then
            log_message "Moving FASTQ files from subdirectory to main directory"
            mv "$srr_accession/fastq"/*.fastq.gz "$srr_accession/" 2>/dev/null || true
        else
            log_message "Error: No FASTQ files found after download for $srr_accession"
            return 1
        fi
    fi
    
    # Check for technical reads (index files) for 10x data
    # Count FASTQ files to determine if we probably have a complete set
    local fastq_count=$(find "$srr_accession" -name "*.fastq.gz" | wc -l)
    
    # For 10x data, we often expect 3 files (v2, with I1) or 4 files (v3, with I1+I2)
    # If we only have 2 files, we might be missing technical reads
    if [ "$fastq_count" -eq 2 ] && [[ "$library_name" == *"10x"* || "$library_name" == *"10X"* || "$library_name" == *"Chromium"* ]]; then
        log_message "Warning: Found only 2 FASTQ files for what appears to be 10x data ($srr_accession)"
        log_message "Index files (I1/I2) may be missing. This might affect cellranger analysis if data is multiplexed."
        
        # Attempt to recover technical reads if enabled
        if [ "$RECOVER_TECHNICAL_READS" = true ]; then
            log_message "Attempting to recover technical reads using fastq-dump..."
            
            # Create temporary directory
            local temp_dir="${srr_accession}/temp_technical"
            mkdir -p "$temp_dir"
            
            # Try using fastq-dump with --include-technical
            if command -v fastq-dump &> /dev/null; then
                log_message "Using fastq-dump to extract technical reads..."
                if ! prefetch "$srr_accession" -O "$temp_dir" &> /dev/null; then
                    log_message "Warning: prefetch failed, continuing with extraction attempt..."
                fi
                
                fastq-dump --split-files --include-technical --gzip -O "$temp_dir" "$srr_accession" &> /dev/null
                
                # Check if we got additional files
                local tech_count=$(find "$temp_dir" -name "*.fastq.gz" | wc -l)
                if [ "$tech_count" -gt "$fastq_count" ]; then
                    log_message "Successfully recovered additional files with fastq-dump"
                    
                    # Move only files that don't already exist
                    for file in "$temp_dir"/*.fastq.gz; do
                        local basename=$(basename "$file")
                        if [ ! -f "$srr_accession/$basename" ]; then
                            mv "$file" "$srr_accession/"
                        fi
                    done
                else
                    log_message "No additional files recovered with fastq-dump"
                fi
                
                # Clean up
                rm -rf "$temp_dir"
            else
                log_message "Warning: fastq-dump not available, skipping technical read recovery"
                rm -rf "$temp_dir"
            fi
        else
            log_message "Continuing with available files. Cellranger should work for single-sample data."
        fi
    fi
    
    # Verify FASTQ integrity
    log_message "Verifying FASTQ integrity for $srr_accession..."
    for fastq in "$srr_accession"/*.fastq.gz; do
        if [ -f "$fastq" ]; then
            # First check if the gzip file is valid
            if ! gzip -t "$fastq" 2>/dev/null; then
                log_message "Warning: FASTQ file $fastq appears to be corrupted"
                continue
            fi
            
            # Validate FASTQ format (check first 1000 entries)
            if ! validate_fastq "$fastq" 1000; then
                log_message "Warning: FASTQ file $fastq has format issues"
                
                # Attempt to fix the FASTQ format if enabled
                if [ "$FIX_FASTQ_FORMAT" = true ]; then
                    if fix_fastq_format "$fastq"; then
                        log_message "Successfully fixed format issues in $fastq"
                    else
                        log_message "Failed to fix format issues in $fastq, file may cause problems with cellranger"
                    fi
                else
                    log_message "FASTQ fixing is disabled. File may cause problems with cellranger."
                fi
            fi
        fi
    done
    
    # Standardize file names
    standardize_fastq_names "$srr_accession"
    
    # Mark as processed
    touch "${STATUS_DIR}/${srr_accession}.processed"
    log_message "Finished processing $srr_accession"
    
    return 0
}

# Function to standardize FASTQ file names
standardize_fastq_names() {
    local srr_accession="$1"
    local file_count=$(find "$srr_accession" -maxdepth 1 -name "*.fastq.gz" | wc -l)
    
    log_message "Standardizing names for $file_count FASTQ files in $srr_accession"
    
    # Get all files sorted
    local files=($(find "$srr_accession" -maxdepth 1 -name "*.fastq.gz" | sort))
    
    # Rename based on file count
    case $file_count in
        1)  # Single-end
            mv "${files[0]}" "$srr_accession/${srr_accession}_1.fastq.gz"
            ;;
        2)  # Paired-end
            mv "${files[0]}" "$srr_accession/${srr_accession}_1.fastq.gz"
            mv "${files[1]}" "$srr_accession/${srr_accession}_2.fastq.gz"
            ;;
        3)  # 10X v2 (I1, R1, R2)
            # Look for index files first
            local i1_file=""
            local r1_file=""
            local r2_file=""
            
            for file in "${files[@]}"; do
                if [[ "$(basename "$file")" == *"_I1_"* ]] || [[ "$(basename "$file")" == *"_index"* ]]; then
                    i1_file="$file"
                elif [[ "$(basename "$file")" == *"_R1_"* ]] || [[ "$(basename "$file")" == *"_1"* ]]; then
                    r1_file="$file"
                elif [[ "$(basename "$file")" == *"_R2_"* ]] || [[ "$(basename "$file")" == *"_2"* ]]; then
                    r2_file="$file"
                fi
            done
            
            # If pattern matching failed, just use positional assignment
            if [ -z "$i1_file" ] || [ -z "$r1_file" ] || [ -z "$r2_file" ]; then
                mv "${files[0]}" "$srr_accession/${srr_accession}_1.fastq.gz"
                mv "${files[1]}" "$srr_accession/${srr_accession}_2.fastq.gz"
                mv "${files[2]}" "$srr_accession/${srr_accession}_3.fastq.gz"
            else
                mv "$i1_file" "$srr_accession/${srr_accession}_1.fastq.gz"
                mv "$r1_file" "$srr_accession/${srr_accession}_2.fastq.gz"
                mv "$r2_file" "$srr_accession/${srr_accession}_3.fastq.gz"
            fi
            ;;
        4)  # 10X v3 (I1, I2, R1, R2)
            mv "${files[0]}" "$srr_accession/${srr_accession}_1.fastq.gz"
            mv "${files[1]}" "$srr_accession/${srr_accession}_2.fastq.gz"
            mv "${files[2]}" "$srr_accession/${srr_accession}_3.fastq.gz"
            mv "${files[3]}" "$srr_accession/${srr_accession}_4.fastq.gz"
            ;;
        *)  # Unexpected number of files
            log_message "Warning: Unexpected number of FASTQ files ($file_count) for $srr_accession"
            log_message "Files may need manual renaming"
            ;;
    esac
    
    log_message "After standardization, files in $srr_accession:"
    ls -l "$srr_accession"/*.fastq.gz 2>/dev/null || log_message "No .fastq.gz files found"
}

# Function to rename to 10X format
rename_to_10x_format() {
    local srr_accession="$1"
    local library_name="$2"
    local lane_number="$3"
    
    # Format lane number with leading zeros
    local lane_num=$(printf "L%03d" "$lane_number")
    
    # Get file count
    local file_count=$(find "${srr_accession}" -maxdepth 1 -name "${srr_accession}_*.fastq.gz" | wc -l)
    
    log_message "Renaming $file_count files from $srr_accession to 10X format with lane $lane_num"
    
    # Create output directory
    mkdir -p "${OUTPUT_DIR}/${library_name}"
    
    # Rename based on file count
    case $file_count in
        4)  # Dual index (I1, I2, R1, R2)
            mv "${srr_accession}/${srr_accession}_1.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_I1_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_2.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_I2_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_3.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R1_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_4.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R2_001.fastq.gz"
            ;;
        3)  # Single index (I1, R1, R2)
            mv "${srr_accession}/${srr_accession}_1.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_I1_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_2.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R1_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_3.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R2_001.fastq.gz"
            ;;
        2)  # Standard paired-end (R1, R2)
            mv "${srr_accession}/${srr_accession}_1.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R1_001.fastq.gz"
            mv "${srr_accession}/${srr_accession}_2.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R2_001.fastq.gz"
            ;;
        1)  # Single-end (R1 only)
            mv "${srr_accession}/${srr_accession}_1.fastq.gz" "${OUTPUT_DIR}/${library_name}/${library_name}_S1_${lane_num}_R1_001.fastq.gz"
            ;;
        *)  # Unexpected number of files
            log_message "Warning: Unexpected number of FASTQ files ($file_count) for $srr_accession"
            return 1
            ;;
    esac
    
    log_message "Files after renaming to 10X format:"
    find "${OUTPUT_DIR}/${library_name}" -name "${library_name}_S1_${lane_num}_*_001.fastq.gz" | sort
    
    return 0
}

# ======== CELLRANGER FUNCTIONS ========
# Function to verify cellranger output
verify_cellranger_output() {
    local sample="$1"
    local cellranger_out="${OUTPUT_DIR}/${sample}/cellranger_outs/outs"
    
    # Check alternate path if standard not found
    if [ ! -d "$cellranger_out" ]; then
        cellranger_out="${OUTPUT_DIR}/${sample}/cellranger_outs/${sample}/outs"
        if [ ! -d "$cellranger_out" ]; then
            log_message "Error: Cellranger output directory not found for $sample"
            return 1
        fi
    fi
    
    # Check for filtered feature-barcode matrix
    if [ ! -d "${cellranger_out}/filtered_feature_bc_matrix" ] || \
       [ ! -f "${cellranger_out}/filtered_feature_bc_matrix/matrix.mtx.gz" ]; then
        log_message "Error: Filtered feature-barcode matrix not found for $sample"
        return 1
    fi
    
    # Check for raw feature-barcode matrix
    if [ ! -d "${cellranger_out}/raw_feature_bc_matrix" ] || \
       [ ! -f "${cellranger_out}/raw_feature_bc_matrix/matrix.mtx.gz" ]; then
        log_message "Error: Raw feature-barcode matrix not found for $sample"
        return 1
    fi
    
    log_message "Cellranger output for $sample verified successfully"
    return 0
}

# Function to clean up cellranger temporary files
clean_cellranger_temp() {
    local sample="$1"
    local sample_dir="${OUTPUT_DIR}/${sample}/cellranger_outs"
    
    # Check if cleanup already done
    if [ -f "${STATUS_DIR}/${sample}.cellranger_cleaned" ]; then
        return 0
    fi
    
    log_message "Cleaning up cellranger temporary files for $sample..."
    
    # Remove large temporary directories if they exist
    for dir in "SC_RNA_COUNTER_CS" "_cmdline" "_filelist" "_finalstate" "_invocation" "_jobmode" "_log" "_mrosource" "_sitecheck" "_tmpdir"; do
        if [ -d "${sample_dir}/${dir}" ]; then
            log_message "Removing ${dir} directory..."
            rm -rf "${sample_dir}/${dir}"
        elif [ -d "${sample_dir}/${sample}/${dir}" ]; then
            log_message "Removing ${dir} from alternate path..."
            rm -rf "${sample_dir}/${sample}/${dir}"
        fi
    done
    
    touch "${STATUS_DIR}/${sample}.cellranger_cleaned"
    log_message "Cleanup completed for $sample"
}

# Function to delete FASTQs after successful cellranger run
delete_sample_fastqs() {
    local sample="$1"
    local sample_dir="${OUTPUT_DIR}/${sample}"
    
    # Skip if keeping FASTQs is requested
    if [ "$KEEP_FASTQS" = true ]; then
        log_message "Keeping FASTQ files for $sample as requested"
        return 0
    fi
    
    # Skip if already deleted
    if [ -f "${STATUS_DIR}/${sample}.fastqs_deleted" ]; then
        return 0
    fi
    
    # Verify cellranger output first
    if verify_cellranger_output "$sample"; then
        log_message "Deleting FASTQ files for $sample after successful cellranger verification..."
        find "$sample_dir" -name "*.fastq.gz" -delete
        touch "${STATUS_DIR}/${sample}.fastqs_deleted"
        log_message "Deleted FASTQ files for $sample"
    else
        log_message "Error: Cellranger output verification failed for $sample, keeping FASTQ files"
        return 1
    fi
}

# Function to run cellranger on a sample
run_cellranger() {
    local sample="$1"
    local sample_dir="${OUTPUT_DIR}/${sample}"
    
    # Skip if already processed
    if [ -f "${STATUS_DIR}/${sample}.cellranger_complete" ]; then
        log_message "Skipping cellranger for $sample - already processed"
        
        # Check if cleanup needed
        if [ ! -f "${STATUS_DIR}/${sample}.cellranger_cleaned" ]; then
            clean_cellranger_temp "$sample"
        fi
        
        # Check if deletion needed
        if [ ! -f "${STATUS_DIR}/${sample}.fastqs_deleted" ] && [ "$KEEP_FASTQS" = false ]; then
            delete_sample_fastqs "$sample"
        fi
        
        return 0
    fi
    
    # Check for required files
    if [ ! -d "$sample_dir" ]; then
        log_message "Error: Sample directory $sample_dir doesn't exist"
        return 1
    fi
    
    local fastq_count=$(find "$sample_dir" -name "*.fastq.gz" | wc -l)
    if [ "$fastq_count" -eq 0 ]; then
        log_message "Error: No FASTQ files found in $sample_dir"
        return 1
    fi
    
    # Validate all FASTQ files before running cellranger
    log_message "Validating all FASTQ files before running cellranger..."
    local validation_failed=false
    
    for fastq in "$sample_dir"/*.fastq.gz; do
        if [ -f "$fastq" ]; then
            if ! validate_fastq "$fastq" 1000; then
                log_message "Warning: FASTQ file $fastq has format issues"
                
                if [ "$FIX_FASTQ_FORMAT" = true ]; then
                    if fix_fastq_format "$fastq"; then
                        log_message "Successfully fixed format issues in $fastq"
                    else
                        log_message "Failed to fix format issues in $fastq"
                        validation_failed=true
                    fi
                else
                    log_message "FASTQ fixing is disabled"
                    validation_failed=true
                fi
            fi
        fi
    done
    
    if [ "$validation_failed" = true ]; then
        log_message "Some FASTQ files failed validation. Continuing with caution."
    fi
    
    # Check for presence of expected file types for 10x data
    local r1_count=$(find "$sample_dir" -name "*_R1_*.fastq.gz" | wc -l)
    local r2_count=$(find "$sample_dir" -name "*_R2_*.fastq.gz" | wc -l)
    local i1_count=$(find "$sample_dir" -name "*_I1_*.fastq.gz" | wc -l)
    local i2_count=$(find "$sample_dir" -name "*_I2_*.fastq.gz" | wc -l)
    
    log_message "Found file counts - R1: $r1_count, R2: $r2_count, I1: $i1_count, I2: $i2_count"
    
    # Warn if we might have incomplete 10x data
    if [[ "$sample" == *"10x"* || "$sample" == *"10X"* || "$sample" == *"Chromium"* ]]; then
        if [ "$r1_count" -eq 0 ] || [ "$r2_count" -eq 0 ]; then
            log_message "Error: Missing required R1 or R2 files for 10x sample $sample. Cellranger requires these files."
            return 1
        fi
        
        if [ "$i1_count" -eq 0 ]; then
            log_message "Warning: No I1 index files found for 10x sample $sample."
            log_message "This might be fine for single-sample data, but could cause issues with multiplexed samples."
        fi
    fi
    
    log_message "Running cellranger for sample $sample with $fastq_count FASTQ files..."
    log_message "FASTQ files: $(find "$sample_dir" -name "*.fastq.gz" | sort)"
    
    # Create output directory
    mkdir -p "${sample_dir}/cellranger_outs"
    
    # Run cellranger with timeout
    if timeout 86400 cellranger count \
        --id="$sample" \
        --transcriptome="$CELLRANGER_REF" \
        --fastqs="$sample_dir" \
        --sample="$sample" \
        --create-bam=false \
        --nosecondary \
        --localcores="$THREADS_PER_JOB" \
        --localmem=128 \
        --output-dir="${sample_dir}/cellranger_outs"; then
        
        # Mark as complete
        touch "${STATUS_DIR}/${sample}.cellranger_complete"
        log_message "Cellranger completed successfully for $sample"
        
        # Verify, clean up, and delete FASTQs
        if verify_cellranger_output "$sample"; then
            clean_cellranger_temp "$sample"
            if [ "$KEEP_FASTQS" = false ]; then
                delete_sample_fastqs "$sample"
            fi
        else
            log_message "Warning: Cellranger output verification failed for $sample"
        fi
    else
        log_message "Error running cellranger for $sample"
        return 1
    fi
}

# ======== MAIN WORKFLOW ========
main() {
    log_message "Starting SRA to cellranger workflow (v${VERSION})"
    
    # Check dependencies and inputs
    check_dependencies
    check_inputs
    
    # Create status directory
    mkdir -p "$STATUS_DIR"
    
    # Display configuration
    log_message "Using $PARALLEL_JOBS parallel jobs with $THREADS_PER_JOB threads each"
    log_message "Output directory: ${OUTPUT_DIR}"
    log_message "Cellranger reference: ${CELLRANGER_REF}"
    log_message "Keep FASTQ files: ${KEEP_FASTQS}"
    log_message "Recover technical reads: ${RECOVER_TECHNICAL_READS}"
    log_message "Fix malformed FASTQ files: ${FIX_FASTQ_FORMAT}"
    
    # Display SRR accessions
    log_message "SRR accessions found in SraRunTable.txt:"
    grep -v "^Run" SraRunTable.txt | cut -f1 | sort | uniq | head -5
    log_message "... and more (total: $(grep -v "^Run" SraRunTable.txt | cut -f1 | sort | uniq | wc -l))"
    
    # Declare associative arrays
    declare -A sample_runs
    declare -A sample_run_counts
    declare -A processed_runs_per_sample
    
    # Step 1: Collect samples and their runs
    log_message "Collecting samples and runs from SraRunTable.txt..."
    
    while IFS=$'\t' read -r run assay_type avgspotlen bases bioproject biosample bytes center consent datastore provider region exp instrument library_name layout selection source organism platform release create ver sample rest || [ -n "$run" ]; do
        if [[ $run != "Run" ]] && [[ -n "$run" ]]; then  # Skip header and empty lines
            if [[ -n "$library_name" ]]; then
                # Store run for this sample
                if [[ -z "${sample_runs[$library_name]}" ]]; then
                    sample_runs[$library_name]="$run"
                    sample_run_counts[$library_name]=1
                else
                    sample_runs[$library_name]="${sample_runs[$library_name]} $run"
                    sample_run_counts[$library_name]=$((sample_run_counts[$library_name] + 1))
                fi
                
                # Initialize processed counter
                if [[ -z "${processed_runs_per_sample[$library_name]}" ]]; then
                    processed_runs_per_sample[$library_name]=0
                fi
            else
                log_message "Warning: Empty library_name for run $run, skipping"
            fi
        fi
    done < SraRunTable.txt
    
    log_message "Found ${#sample_runs[@]} unique samples"
    
    # Step 2: Process SRA runs in parallel
    if command -v parallel &> /dev/null; then
        log_message "Using GNU parallel to process SRA runs (max $PARALLEL_JOBS jobs)..."
        
        # Create temporary file with runs to process
        TMP_RUNLIST=$(mktemp)
        while IFS=$'\t' read -r run assay_type avgspotlen bases bioproject biosample bytes center consent datastore provider region exp instrument library_name layout selection source organism platform release create ver sample rest || [ -n "$run" ]; do
            if [[ $run != "Run" ]] && [[ -n "$run" ]] && [[ -n "$library_name" ]]; then
                echo "$run $library_name" >> "$TMP_RUNLIST"
            fi
        done < SraRunTable.txt
        
        # Export functions and variables for parallel
        export -f process_sra_dir check_fastq_files_exist standardize_fastq_names log_message validate_fastq fix_fastq_format
        export THREADS_PER_JOB STATUS_DIR OUTPUT_DIR DOWNLOAD_TIMEOUT FIX_FASTQ_FORMAT RECOVER_TECHNICAL_READS
        
        parallel -j "$PARALLEL_JOBS" --colsep ' ' 'process_sra_dir {1} {2}' :::: "$TMP_RUNLIST"
        rm -f "$TMP_RUNLIST"
    else
        log_message "GNU parallel not found, processing SRA runs sequentially..."
        while IFS=$'\t' read -r run assay_type avgspotlen bases bioproject biosample bytes center consent datastore provider region exp instrument library_name layout selection source organism platform release create ver sample rest || [ -n "$run" ]; do
            if [[ $run != "Run" ]] && [[ -n "$run" ]] && [[ -n "$library_name" ]]; then
                process_sra_dir "$run" "$library_name" 
            fi
        done < SraRunTable.txt
    fi
    
    # Step 3: Organize files by sample and run cellranger
    log_message "Organizing files by sample and running cellranger..."
    
    for sample in "${!sample_run_counts[@]}"; do
        log_message "Processing sample: $sample"
        
        # Skip if cellranger already complete
        if [ -f "${STATUS_DIR}/${sample}.cellranger_complete" ]; then
            log_message "Cellranger already complete for $sample"
            
            # Check if cleanup or deletion needed
            if [ ! -f "${STATUS_DIR}/${sample}.cellranger_cleaned" ]; then
                clean_cellranger_temp "$sample"
            fi
            
            if [ ! -f "${STATUS_DIR}/${sample}.fastqs_deleted" ] && [ "$KEEP_FASTQS" = false ]; then
                delete_sample_fastqs "$sample"
            fi
            
            continue
        fi
        
        # Create sample directory
        mkdir -p "${OUTPUT_DIR}/${sample}"
        
        # Organize files if not already done
        if [ ! -f "${STATUS_DIR}/${sample}.10x_formatted" ]; then
            log_message "Converting FASTQs for sample $sample to 10X format..."
            
            # Get all runs for this sample
            read -r -a runs <<< "${sample_runs[$sample]}"
            # Do not use "local" here as we're not in a function
            lane=1
            
            for run in "${runs[@]}"; do
                if [ -f "${STATUS_DIR}/${run}.processed" ]; then
                    log_message "Organizing run $run as lane $lane for sample $sample"
                    
                    if rename_to_10x_format "$run" "$sample" "$lane"; then
                        processed_runs_per_sample[$sample]=$((processed_runs_per_sample[$sample] + 1))
                    else
                        log_message "Error organizing run $run for sample $sample"
                    fi
                    
                    lane=$((lane + 1))
                else
                    log_message "Skipping run $run for sample $sample - not yet processed"
                fi
            done
            
            # Mark as 10X formatted if all runs were processed
            if [ "${processed_runs_per_sample[$sample]}" -eq "${sample_run_counts[$sample]}" ]; then
                touch "${STATUS_DIR}/${sample}.10x_formatted"
                log_message "All runs for $sample have been converted to 10X format"
            fi
        fi
        
        # Run cellranger if all runs are processed
        if [ "${processed_runs_per_sample[$sample]}" -eq "${sample_run_counts[$sample]}" ]; then
            log_message "All runs for sample $sample have been processed and organized"
            run_cellranger "$sample"
        else
            log_message "Still waiting for $((sample_run_counts[$sample] - processed_runs_per_sample[$sample])) runs to be processed for sample $sample"
        fi
    done
    
    # Clean up empty directories
    log_message "Cleaning up empty SRR directories..."
    find "$BASE_DIR" -type d -name "SRR*" -empty -delete
    
    # Print final status report
    log_message "\nProcessing Status:"
    log_message "===================="
    log_message "SRA Files processed: $(find "${STATUS_DIR}" -name "*.processed" | wc -l)"
    log_message "Cellranger Analyses Completed: $(find "${STATUS_DIR}" -name "*.cellranger_complete" | wc -l)"
    log_message "Cellranger Temp Files Cleaned: $(find "${STATUS_DIR}" -name "*.cellranger_cleaned" | wc -l)"
    log_message "FASTQ Files Deleted: $(find "${STATUS_DIR}" -name "*.fastqs_deleted" | wc -l)"
    
    # List successfully processed samples
    log_message "\nProcessed Samples:"
    for f in "${STATUS_DIR}"/*.cellranger_complete; do
        if [ -f "$f" ]; then
            sample=$(basename "$f" .cellranger_complete)
            fastq_status="FASTQs present"
            if [ -f "${STATUS_DIR}/${sample}.fastqs_deleted" ]; then
                fastq_status="FASTQs deleted"
            fi
            log_message "  - $sample ($fastq_status)"
        fi
    done
    
    log_message "\nScript execution completed."
}

# Setup cleanup trap
trap 'log_message "Script interrupted. Cleaning up temporary directories..."; rm -rf "${BASE_DIR}"/*/tmp; log_message "Cleanup complete."' EXIT

# Start the workflow
main "$@"
