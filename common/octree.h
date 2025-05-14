#ifndef OCTREE_H
#define OCTREE_H

#include <cstdint>
#include <string>
#include <vector>

typedef uint64_t Key;

typedef struct {
  float x, y, z, vx, vy, vz, m; // each body now has a mass 'm'
  Key key;
  int index; // index of the body in the bodies array
} Body;

Key computeMortonKey(Body b, int resolution);
int binaryLength(Key n);

template <typename T> std::string binaryString(T val) {
  std::string s = "";
  if (val == 0)
    return "0";
  while (val > 0) {
    s = (val % 2 == 0 ? "0" : "1") + s;
    val /= 2;
  }
  return s;
}

class Node {
public:
  Key key;
  Node *children[8];
  uint8_t maskChildren;
  std::vector<int> bodyIdx; // indeces into the bodies array
  uint64_t subTreeSize;     // used to build autoropes

  float cx, cy, cz;         // centroid pos
  float tm;                 // total mass

  Node(Key key);
  ~Node();
};

typedef struct {
  int index;
  Key key;
  bool isLeaf;
  int autorope; // index of the next node in the autorope

  float x, y, z; // center of mass of the node
  float mass;    // total mass of the node
  int nBodies;

  std::vector<int> bodies; // indeces into the bodies array
} DFTNode;

class Octree {

public:
  int resolution;
  Node *root;

  Octree(int resolution);
  ~Octree();
  bool isEmpty(Node *&n);
  bool isLeaf(Node *&n);
  int getOctantIndex(Key key, int level);
  void updateAggregateStats(Node *&node, Body &body);
  int subdivide(Node *&node, int level);
  void insert(Body &body);
  void printTree(Node *&node, int level);
  void buildDFT(std::vector<DFTNode> &nodes, Body *bodies);
  void traverse(Node *node, std::vector<DFTNode> &nodes);
};

#endif // OCTREE_H
