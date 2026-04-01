#!/usr/bin/env python3
"""
Full signal extraction script for ninetails package.
Extracts the complete raw signal from a POD5 file for a single read,
along with calibration parameters needed for rescaling to picoamps.

Used by plot_squiggle_pod5() and plot_tail_range_pod5() functions.
"""

import argparse
import pod5
import numpy as np
import pickle
import sys


def winsorize_signal(signal, lower_quantile=0.002, upper_quantile=0.998):
    """Apply winsorization to signal - matching R implementation.
    
    Parameters
    ----------
    signal : numpy.ndarray
        Raw signal array
    lower_quantile : float
        Lower quantile for clipping (default 0.002 = 0.2%)
    upper_quantile : float
        Upper quantile for clipping (default 0.998 = 99.8%)
    
    Returns
    -------
    numpy.ndarray
        Winsorized signal as integers
    """
    if len(signal) == 0:
        return signal
    
    # Match R implementation: use specified percentiles
    lower = np.percentile(signal, lower_quantile * 100)
    upper = np.percentile(signal, upper_quantile * 100)
    
    # Clip values
    signal = np.clip(signal, lower, upper)
    
    # Convert to integer (matching R as.integer)
    return signal.astype(int)


def main():
    parser = argparse.ArgumentParser(
        description="Extract full signal from POD5 file for a single read"
    )
    parser.add_argument(
        "--read_id", 
        required=True, 
        help="Read ID to extract"
    )
    parser.add_argument(
        "--pod5_file", 
        required=True, 
        help="Path to POD5 file"
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Output pickle file for signal and metadata"
    )
    parser.add_argument(
        "--winsorize", 
        action="store_true", 
        default=False, 
        help="Apply winsorization to signal"
    )
    
    args = parser.parse_args()
    
    read_id = str(args.read_id).strip()
    result = None
    
    try:
        with pod5.Reader(args.pod5_file) as reader:
            for read in reader.reads():
                current_id = str(read.read_id).strip()
                
                if current_id == read_id:
                    # Get raw signal
                    signal = read.signal
                    
                    # Get calibration parameters
                    # pA = (signal + offset) * scale
                    calibration_offset = read.calibration.offset
                    calibration_scale = read.calibration.scale
                    
                    # Get sample rate from run_info
                    sample_rate = read.run_info.sample_rate
                    
                    if args.winsorize:
                        signal = winsorize_signal(signal)
                    else:
                        signal = np.array(signal)
                    
                    result = {
                        "read_id": read_id,
                        "signal": signal.tolist(),
                        "signal_length": len(signal),
                        "calibration_offset": float(calibration_offset),
                        "calibration_scale": float(calibration_scale),
                        "sample_rate": int(sample_rate)
                    }
                    
                    break
        
        if result is None:
            print(
                f"Error: Read ID '{read_id}' not found in {args.pod5_file}", 
                file=sys.stderr
            )
            sys.exit(1)
        
        # Save result
        with open(args.output, "wb") as f:
            pickle.dump(result, f)
        
        print(
            f"Extracted signal for {read_id}: {result['signal_length']} samples, "
            f"sample_rate={result['sample_rate']} Hz"
        )
    
    except Exception as e:
        print(f"Error processing POD5 file: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
