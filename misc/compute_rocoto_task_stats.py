"""
Compute Statistics of Task Times From Rocoto

Usage:
    rocotostat -w <WORKFLOW>.xml -d <WORKFLOW>.db > rocoto.out
    python compute_rocoto_task_stats.py rocoto.out

shawn.s.murdzek@noaa.gov
Original draft: Gemini
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import argparse
import sys


#---------------------------------------------------------------------------------------------------
# Program
#---------------------------------------------------------------------------------------------------

def compute_task_stats(filename):
    try:
        # Load the data
        df = pd.read_csv(
            filename,
            sep=r'\s+',
            skiprows=2,
            comment='=',
            names=["CYCLE", "TASK", "JOBID", "STATE", "EXIT_STATUS", "TRIES", "DURATION"],
            na_values=['-']      # Automatically converts '-' strings into missing values (NaN)
        )
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        sys.exit(1)

    # Filter 1: Drop rows where JOBID is -1 and STATE is not 'SUCCEEDED'
    df = df[df['JOBID'] != -1]
    df = df[df['STATE'] == 'SUCCEEDED']

    # Filter 2: Drop rows where DURATION is missing
    df = df.dropna(subset=['DURATION'])

    # Group by TASK, calculate the mean DURATION, and sort from longest to shortest
    task_stats = (
        df.groupby('TASK')['DURATION']
        .agg(['mean', 'std'])
        .round(2)
        .sort_values(by='mean', ascending=False)
    )

    # If a task only appears once, its standard deviation is undefined (NaN).
    # We fill these with 0.0 to keep the output table clean.
    task_stats['std'] = task_stats['std'].fillna(0.0)

    # Print the formatted results
    print(f"{'TASK':<30} | {'AVERAGE DURATION':<15} | {'STD DEVIATION'}")
    print("-" * 65)

    if task_stats.empty:
        print("No valid tasks found after filtering.")
    else:
        # Iterate through the DataFrame rows to print cleanly
        for task, row in task_stats.iterrows():
            print(f"{task:<30} | {row['mean']:<15.2f} | {row['std']:.2f}")

if __name__ == "__main__":
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(description="Compute task statistics from a Rocoto log file.")
    parser.add_argument("filename", help="The path to the text file containing the log output.")

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Run the function with the provided filename
    compute_task_stats(args.filename)


"""
End compute_rocoto_task_stats.py
"""
