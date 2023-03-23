
import argparse

parser = argparse.ArgumentParser(description='python script to test argument parser')

parser.add_argument('-d', default=5, help='some random parameter')

print(
