#!/usr/bin/env python3
"""
Simple POD5 signal extraction script for ninetails package.
This script extracts only the raw signals from POD5 files.
Pseudomove computation is handled on the R side.

In this script, only the polyA tail signal is extracted
based on coordinates provided in corresponding 
Dorado summary file.
"""

import argparse
import pandas as pd
import pod5
import numpy as np
from multiprocessing import Pool
import pickle
import os
import sys

def winsorize_signal(signal):
    """Apply winsorization to signal - matching R implementation"""
    if len(signal) == 0:
        return signal
    
    # Match R implementation: use 0.5% and 99.5% percentiles
    # see ?ninetails::winsorize_signal() in R
    lower = np.percentile(signal, 0.5)  # 0.005 * 100
    upper = np.percentile(signal, 99.5)  # 0.995 * 100
    
    # Clip values
    signal = np.clip(signal, lower, upper)
    
    # Convert to integer (matching R as.integer)
    return signal.astype(int)

def process_pod5_file(args):
    """Process a single POD5 file and extract signals"""
    filename, reads_data, pod5_dir = args
    pod5_path = os.path.join(pod5_dir, filename)
    
    if not os.path.exists(pod5_path):
        print(f"Warning: POD5 file not found: {pod5_path}", file=sys.stderr)
        return {}
    
    signals = {}
    
    try:
        with pod5.Reader(pod5_path) as reader:
            # Create read_id to read mapping (convert UUID to string)
            read_dict = {}
            for read in reader.reads():
                read_id_str = str(read.read_id).strip()
                read_dict[read_id_str] = read
            
            for _, row in reads_data.iterrows():
                read_id = str(row["read_id"]).strip()
                if read_id not in read_dict:
                    continue
                
                read = read_dict[read_id]
                signal = read.signal
                
                start_idx = int(row["poly_tail_start"])
                end_idx = int(row["poly_tail_end"])
                
                if (start_idx > 0 and end_idx > 0 and 
                    start_idx < end_idx and end_idx <= len(signal)):
                    
                    polya_signal = signal[start_idx:end_idx]
                    
                    # First winsorize (this already returns integers)
                    polya_signal = winsorize_signal(polya_signal)
                    
                    # Then interpolate (matching R implementation)
                    n_points = int(np.ceil(0.2 * len(polya_signal)))
                    if n_points > 0:
                        x = np.linspace(0, len(polya_signal) - 1, n_points)
                        polya_signal = np.interp(x, 
                                               np.arange(len(polya_signal)), 
                                               polya_signal)
                        # Round and convert to int to match R behavior
                        polya_signal = np.round(polya_signal).astype(int)
                    
                    signals[read_id] = polya_signal.tolist()
                else:
                    signals[read_id] = []
    
    except Exception as e:
        print(f"Error processing {pod5_path}: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
    
    return signals

def main():
    parser = argparse.ArgumentParser(description="Extract polyA signals from POD5 files")
    parser.add_argument("--input", required=True, help="Input CSV file with read information")
    parser.add_argument("--pod5_dir", required=True, help="Directory containing POD5 files")
    parser.add_argument("--output", required=True, help="Output pickle file for signals")
    parser.add_argument("--num_cores", type=int, default=1, help="Number of CPU cores to use")
    
    args = parser.parse_args()
    
    # Read input data
    print(f"Reading input data from {args.input}")
    polya_data = pd.read_csv(args.input)
    print(f"Found {len(polya_data)} reads to process")
    
    # Ensure read_id is string type
    polya_data["read_id"] = polya_data["read_id"].astype(str)
    
    # Group by filename
    grouped = polya_data.groupby("filename")
    file_count = len(grouped)
    print(f"Processing {file_count} POD5 files with {args.num_cores} cores")
    
    # Prepare arguments for parallel processing
    process_args = [(filename, group, args.pod5_dir) 
                    for filename, group in grouped]
    
    # Process in parallel or sequentially
    if args.num_cores > 1:
        with Pool(args.num_cores) as pool:
            results = pool.map(process_pod5_file, process_args)
    else:
        results = [process_pod5_file(arg) for arg in process_args]
    
    # Combine results
    all_signals = {}
    for result in results:
        all_signals.update(result)
    
    print(f"Extracted signals for {len(all_signals)} reads")
    
    # Save results
    with open(args.output, "wb") as f:
        pickle.dump(all_signals, f)
    
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()
