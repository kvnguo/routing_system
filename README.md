# Shortest Path Finder for Cities
Author: Kevin Guo  
This program finds the shortest path between two cities based on distance or time.

## Description
This C++ program computes and visualizes the shortest path between cities stored in a dataset. It utilizes Dijkstra's algorithm to find the path based on either geographical distance or travel time, depending on the user's choice. The program reads city data from a file, calculates distances between cities using geographic coordinates, and categorizes cities into different types (LOCAL, METRO, REGIONAL, NATIONAL) based on population and proximity.

## Features
- **Input**: Reads city data from a text file (`cities.txt`) including city name, type, latitude, longitude, and population.
- **Output**: Generates tables and files (`vertex_table.txt`, `dist_table.txt`, `time_table.txt`, `edge_table.txt`) with detailed information about cities, distances, and connections.
- **Functionality**: 
  - Computes distances between all cities based on geographic coordinates.
  - Categorizes cities into types (LOCAL, METRO, REGIONAL, NATIONAL) and updates city types based on population and distance criteria.
  - Finds the shortest path between two cities using Dijkstra's algorithm, considering either geographical distance or travel time.

## Usage
### Compilation
g++ -std=c++11 RoutingSystem.cpp -o RoutingSystem

### Execution
1. Run the application with the desired mode and input file:
   ```
   ./RoutingSystem -info|dist|time cities.csv
   ```
    - `-info`: Generates detailed tables of cities, distances, times, and connections.
    - `-dist`: Computes shortest path based on geographical distance.
    - `-time`: Computes shortest path based on travel time (by truck for non-NATIONAL cities, and by plane for NATIONAL cities).

2. In routing mode, enter the source and destination cities when prompted:
   ```
   Enter> [source_city] [destination_city]
   ```
   - `Ctrl+D` or `Ctrl+Z` (Windows) to exit the application

## Files Generated
- **vertex_table.txt**: Detailed information about each city.
- **dist_table.txt**: Distance matrix between cities.
- **time_table.txt**: Travel time matrix between cities.
- **edge_table.txt**: Connectivity details showing adjacent cities for each city.




