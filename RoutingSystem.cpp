/*
Kevin Guo
Lab 6: This program finds the shortest path between two cities based on distance or time.
*/

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stack>
#include <unordered_map>
#include <set>

using namespace std;

// Constants and enums
#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

const float earthradius = 3963.0;    // [miles]
const float distance_accuracy = 5.0; // [miles]

const int national_minpop = 1000000;

const float national_dist = 150.0; // [miles]
const float regional_dist = 100.0; // [miles]

const float local_maxdist = 50.0; // [miles]

const float plane_speed = 520.0; // [mph]
const float truck_speed = 65.0;  // [mph]

enum city_t
{
  LOCAL,
  METRO,
  REGIONAL,
  NATIONAL
};

enum color_t
{
  WHITE,
  BLACK
};

// Represents a city with its name, type, latitude, longitude, and population
class city
{
public:
  static int nameWidth;

  // < overload for sorting by population
  bool operator<(const city &rhs) const { return rhs.population < population; }

  // Getters
  string get_name() const { return name; }
  city_t get_type() const { return type; }
  float get_latitude() const { return latitude; }
  float get_longitude() const { return longitude; }
  int get_population() const { return population; }

  // Setters
  void set_name(const string name) { this->name = name; }
  void set_type(const city_t type) { this->type = type; }
  void set_latitude(const float latitude) { this->latitude = latitude; }
  void set_longitude(const float longitude) { this->longitude = longitude; }
  void set_population(const int population) { this->population = population; }

private:
  string name;
  city_t type;
  float latitude, longitude;
  int population;
};
int city::nameWidth = 0;

// Takes in and sets city data
istream &operator>>(istream &in, city &c)
{
  static unordered_map<string, city_t> string2city_t = {{"LOCAL", LOCAL}, {"METRO", METRO}, {"REGIONAL", REGIONAL}, {"NATIONAL", NATIONAL}};

  string name, type;
  float latitude, longitude;
  int population;
  in >> name >> type >> latitude >> longitude >> population;

  c.set_name(name);
  city::nameWidth = max(static_cast<size_t>(city::nameWidth), name.size() + 3);
  c.set_type(string2city_t[type]);
  c.set_latitude(static_cast<float>(latitude * DEG2RAD));
  c.set_longitude(static_cast<float>(longitude * DEG2RAD));
  c.set_population(population);

  return in;
}

// Prints out city data
ostream &operator<<(ostream &out, const city &c)
{
  static unordered_map<city_t, string> city_t2string = {{LOCAL, "LOCAL"}, {METRO, "METRO"}, {REGIONAL, "REGIONAL"}, {NATIONAL, "NATIONAL"}};

  out << setw(city::nameWidth) << left << setfill('.') << c.get_name() << setfill(' ')
      << "  " << setw(8) << left << city_t2string[c.get_type()]
      << "  " << setw(8) << right << c.get_population()
      << "  " << setw(7) << fixed << setprecision(2) << right << c.get_latitude() * RAD2DEG
      << "  " << setw(7) << fixed << setprecision(2) << right << c.get_longitude() * RAD2DEG;
  return out;
}

// Represents a symmetric matrix using a vector
class matrix
{
public:
  matrix(int n) : data(n * (n + 1) / 2) {}

  // () overload to access and store data
  float &operator()(int x, int y)
  {
    // For symmetry
    if (y > x)
      swap(x, y);

    // Row calculations + column offset
    return data[x * (x + 1) / 2 + y];
  }

private:
  vector<float> data;
};

// Creates a vertex table from a list of city
void create_vertex_table(string fname, vector<city> &vertex_table)
{
  ifstream fin(fname);
  string line;
  while (getline(fin, line))
  {
    // Correct formmatting
    replace(line.begin(), line.end(), ',', ' ');
    auto first_comma = std::find(line.begin(), line.end(), ' ');
    *first_comma = '_';

    // Create and insert city
    city newCity;
    istringstream sin(line);
    sin >> newCity;
    vertex_table.push_back(newCity);
  }
  std::sort(vertex_table.begin(), vertex_table.end()); // Sort cities by populations
  fin.close();
}

// Updates cites in vertex table based on label
void update_vertex_table(vector<city> &vertex_table, matrix &dist_table)
{
  // Relabel METRO cities as NATIONAL or REGIONAL
  for (auto &city : vertex_table)
  {
    if (city.get_type() == METRO)
      city.set_type(city.get_population() > national_minpop ? NATIONAL : REGIONAL);
  }

  // Only largest NATIONAL city within 150 miles gets to keep the title
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    int row = i - 1;
    for (size_t j = 1; j <= i; j++)
    {
      int col = j - 1;

      if (vertex_table[row].get_type() == NATIONAL && vertex_table[col].get_type() == NATIONAL)
        if (dist_table(row, col) < national_dist)
          if (vertex_table[row].get_population() < vertex_table[col].get_population())
            vertex_table[row].set_type(REGIONAL);
    }
  }

  // Only largest REGIONAL city within 100 miles gets to keep the title
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    int row = i - 1;
    for (size_t j = 1; j <= i; j++)
    {
      int col = j - 1;
      if (vertex_table[row].get_type() == REGIONAL && vertex_table[col].get_type() == REGIONAL)
        if (dist_table(row, col) < regional_dist)
          if (vertex_table[row].get_population() < vertex_table[col].get_population())
            vertex_table[row].set_type(LOCAL);
    }
  }
}

// Creates table that stores distances between all cities
void create_dist_table(vector<city> &vertex_table, matrix &dist_table)
{
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    int row = i - 1;
    city city_a = vertex_table[row];
    float lat_a = city_a.get_latitude();

    for (size_t j = 1; j <= i; j++)
    {
      int col = j - 1;

      city city_b = vertex_table[col];

      float lat_diff = city_b.get_latitude() - city_a.get_latitude();
      float lon_diff = city_b.get_longitude() - city_a.get_longitude();

      float lat_b = city_b.get_latitude();

      float central_angle = pow(sin(lat_diff / 2), 2) + cos(lat_a) * cos(lat_b) * pow(sin(lon_diff / 2), 2);

      float angular_distance = 2 * asin(sqrt(central_angle));

      float distance = earthradius * angular_distance;
      distance = round(distance / distance_accuracy) * distance_accuracy;

      dist_table(row, col) = distance;
    }
  }
}

// Creates table that stores time between all cities
void create_time_table(vector<city> &vertex_table, matrix &dist_table, matrix &time_table)
{
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    int row = i - 1;
    for (size_t j = 1; j <= i; j++)
    {
      int col = j - 1;

      // NATIONAL city tranversal is done by plane
      if (vertex_table[row].get_type() == NATIONAL && vertex_table[col].get_type() == NATIONAL)
        time_table(row, col) = dist_table(row, col) / plane_speed;
      // All other travel is done by truck
      else
        time_table(row, col) = dist_table(row, col) / truck_speed;
    }
  }
}

// Simple representation of a nonlocal city with its index and distance
struct nonlocal
{
  int index;
  float distance;

  nonlocal(int index, float distance) : index(index), distance(distance) {}

  // < overload for sorting by distance
  bool operator<(const nonlocal &rhs) const { return distance < rhs.distance; }
};

// Creates an adjacency matrix
void create_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table)
{
  for (size_t i = 0; i < vertex_table.size(); i++)
  {
    // All NATIONAL cities connect to each other
    if (vertex_table[i].get_type() == NATIONAL)
    {
      for (size_t j = i + 1; j < vertex_table.size(); j++)
      {
        if (vertex_table[j].get_type() == NATIONAL)
        {
          edge_table[i].insert(j);
          edge_table[j].insert(i);
        }
      }
    }
    // REGIONAL cities connect to the 3 closest NONLOCAL cities
    else if (vertex_table[i].get_type() == REGIONAL)
    {
      vector<nonlocal> distances;
      for (size_t j = 0; j < vertex_table.size(); j++)
      {
        if (i != j && vertex_table[j].get_type() != LOCAL)
          distances.emplace_back(j, dist_table(i, j));
      }

      partial_sort(distances.begin(), distances.begin() + 3, distances.end());

      for (size_t j = 0; j < static_cast<size_t>(3); j++)
      {
        edge_table[i].insert(distances[j].index);
        edge_table[distances[j].index].insert(i);
      }
    }
    // Local cities connect to the 5 closest NONLOCAL cities and local cities within 50 miles
    else
    {
      vector<nonlocal> distances;
      for (size_t j = 0; j < vertex_table.size(); j++)
      {
        if (i != j)
        {
          if (vertex_table[j].get_type() != LOCAL)
            distances.emplace_back(j, dist_table(i, j));
          else if (dist_table(i, j) < local_maxdist)
          {
            edge_table[i].insert(j);
            edge_table[j].insert(i);
          }
        }
      }

      partial_sort(distances.begin(), distances.begin() + 5, distances.end());

      for (size_t j = 0; j < static_cast<size_t>(5); j++)
      {
        edge_table[i].insert(distances[j].index);
        edge_table[distances[j].index].insert(i);
      }
    }
  }
}

// Writes vertex_table to file
void write_vertex_table(vector<city> &vertex_table)
{
  ofstream fout("vertex_table.txt");
  ostringstream oss;
  int index = 0;
  for (auto &city : vertex_table)
    oss << right << setw(4) << index++ << "  " << city << endl;
  fout << oss.str();
  fout.close();
}

// Writes dist_table to file
void write_dist_table(vector<city> &vertex_table, matrix &dist_table)
{
  ofstream fout("dist_table.txt");
  ostringstream oss;
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    int row = i - 1;
    for (size_t j = 1; j <= i; j++)
    {
      int col = j - 1;

      if (row != col)
      {
        oss << right << setw(4) << row << "  " << setw(city::nameWidth) << left
            << setfill('.') << vertex_table[row].get_name() << " to " << setw(city::nameWidth)
            << vertex_table[col].get_name() << right << setfill(' ') << setw(8)
            << fixed << setprecision(1) << dist_table(row, col) << " miles" << endl;
      }
    }
    if (i != 1)
      oss << endl;
  }
  fout << oss.str();
  fout.close();
}

// Writes time_table to file
void write_time_table(vector<city> &vertex_table, matrix &time_table)
{
  ofstream fout("time_table.txt");
  ostringstream oss;
  for (size_t i = 1; i <= vertex_table.size(); i++)
  {
    for (size_t j = 1; j <= i; j++)
    {
      int row = i - 1;
      int col = j - 1;

      if (row != col)
      {
        oss << right << setw(4) << row << "  " << setw(city::nameWidth) << left
            << setfill('.') << vertex_table[row].get_name() << " to " << setw(city::nameWidth)
            << vertex_table[col].get_name() << right << setfill(' ') << setw(8)
            << fixed << setprecision(1) << time_table(row, col) << " hours" << endl;
      }
    }
    if (i != 1)
      oss << endl;
  }
  fout << oss.str();
  fout.close();
}

// Writes edge_table to file
void write_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table, matrix &time_table)
{
  ofstream fout("edge_table.txt");
  ostringstream oss;
  for (size_t row = 0; row < edge_table.size(); row++)
  {
    oss << right << setw(4) << row << " " << vertex_table[row].get_name() << endl;
    for (auto &col : edge_table[row])
    {
      oss << "  " << right << setw(4) << col << ' ' << setw(city::nameWidth) << left
          << setfill('.') << vertex_table[col].get_name() << right << setfill(' ')
          << setw(8) << fixed << setprecision(1) << dist_table(row, col) << " miles"
          << setw(5) << time_table(row, col) << " hours" << endl;
    }
    if (row != edge_table.size() - 1)
      oss << endl;
  }
  fout << oss.str();
  fout.close();
}

// Finds the shortest path between two cities based on distance or time and shows route 
void dijkstra_route(int city_from, int city_to, vector<city> &vertex_table, vector<set<int>> &edge_table, char mode, matrix &dist_table, matrix &time_table)
{
  vector<color_t> vcolor(vertex_table.size(), WHITE);
  vector<float> vdist(vertex_table.size(), numeric_limits<float>::max());
  vector<int> vlink(vertex_table.size(), -1);

  vdist[city_from] = 0; // Set starting city distance to 0

  while (true)
  {
    // Find unvisited vertex with minimum distance
    int next_i = -1;
    float mindist = numeric_limits<float>::max();
    for (size_t i = 0; i < vcolor.size(); i++)
    {
      if (vcolor[i] == WHITE && mindist > vdist[i])
      {
        next_i = i;
        mindist = vdist[i];
      }
    }

    int i = next_i;
    if (i == -1 || i == city_to)
      break; // No unvisited vertices or destination reached

    vcolor[i] = BLACK; // Mark current vertex as visited

    // Update distances and links for neighboring vertices
    for (auto j : edge_table[i])
    {
      float wij = mode == 'd' ? dist_table(i, j) : time_table(i, j);
      if (vcolor[j] == WHITE && vdist[j] > vdist[i] + wij)
      {
        vdist[j] = vdist[i] + wij;
        vlink[j] = i;
      }
    }
  }

  if (vdist[city_to] == numeric_limits<float>::max())
  {
    cout << "No path found\n"; // No path exists
  }
  else
  {
    stack<int> S;
    // Build path from destination to starting city
    for (int i = city_to; i != city_from; i = vlink[i])
      S.push(i);
    S.push(city_from);

    int prev = S.top();
    float total_dist = 0;
    float total_time = 0;

    // Print starting city information
    cout << right << setw(8) << fixed << setprecision(2) << total_dist << " miles "
         << setw(5) << fixed << total_time << " hours  "
         << setw(3) << prev << ' '
         << setw(city::nameWidth) << left << setfill('.') << vertex_table[prev].get_name()
         << right << setw(10) << setfill(' ') << vertex_table[prev].get_population() << endl;

    S.pop();

    // Print path information
    while (!S.empty())
    {
      int curr = S.top();
      total_dist += dist_table(prev, curr);
      total_time += time_table(prev, curr);

      cout << right << setw(8) << fixed << setprecision(2) << total_dist << " miles "
           << setw(5) << fixed << total_time << " hours  "
           << setw(3) << curr << ' '
           << setw(city::nameWidth) << left << setfill('.') << vertex_table[curr].get_name()
           << right << setw(10) << setfill(' ') << vertex_table[curr].get_population()
           << setw(8) << fixed << setprecision(2) << dist_table(prev, curr) << " miles  "
           << setw(4) << fixed << time_table(prev, curr) << " hours" << endl;

      prev = S.top();
      S.pop();
    }
    cout << endl;
  }
}

int main(int argc, char *argv[])
{
  // Parse commandline arguments
  bool info = false, dist = false, time = false;
  if (argc == 3)
  {
    string flag(argv[1]);
    info = flag == "-info";
    dist = flag == "-dist";
    time = flag == "-time";

    if (!(info || dist || time))
    {
      cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
      return 1;
    }
  }
  else
  {
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    return 1;
  }

  // Create vertex_table
  vector<city> vertex_table;
  create_vertex_table(argv[argc - 1], vertex_table);

  // Create dist_table
  matrix dist_table(vertex_table.size());
  create_dist_table(vertex_table, dist_table);

  // Update vertex_table
  update_vertex_table(vertex_table, dist_table);

  // Create time_table
  matrix time_table(vertex_table.size());
  create_time_table(vertex_table, dist_table, time_table);

  // Create edge_table
  vector<set<int>> edge_table(vertex_table.size());
  create_edge_table(vertex_table, edge_table, dist_table);

  if (info)
  {
    write_vertex_table(vertex_table);
    write_dist_table(vertex_table, dist_table);
    write_time_table(vertex_table, time_table);
    write_edge_table(vertex_table, edge_table, dist_table, time_table);
  }
  else
  {
    while (true)
    {
      cout << "Enter> ";

      string city_from, city_to;
      int city_from_index = -1, city_to_index = -1;
      cin >> city_from >> city_to;

      if (cin.eof())
        break;

      // Checking is city exist 
      for (size_t i = 0; i < vertex_table.size(); i++)
      {
        if (vertex_table[i].get_name() == city_from)
          city_from_index = i;
        if (vertex_table[i].get_name() == city_to)
          city_to_index = i;
        if (city_from_index != -1 && city_to_index != -1)
          break;
      }

      // Error messages
      if (city_from_index == -1)
      {
        cout << city_from << ": prefix not found\n"
             << endl;
        continue;
      }
      if (city_to_index == -1)
      {
        cout << city_to << ": prefix not found\n"
             << endl;
        continue;
      }

      dijkstra_route(city_from_index, city_to_index, vertex_table, edge_table, dist ? 'd' : 't', dist_table, time_table);
    }
    cout << endl;
  }
}
