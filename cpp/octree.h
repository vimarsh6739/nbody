#ifndef OCTREE_H
#define OCTREE_H

#include <string>
#include <vector>

typedef long long Key;

typedef struct {
  float x, y, z, vx, vy, vz, m; // each body now has a mass 'm'
  Key key;
  int index; // index of the body in the bodies array
} Body;

Key getKeyNoPrepend(Body b);
int binaryLength(Key n);
std::string binaryString(Key key);

class Node {
public:
  Key key;
  bool isLeaf;

  Node *parent;
  Node *children[8];
  char whichChildren;

  std::vector<int> bodies; // indeces into the bodies array

  int subTreeSize; // used to build autoropes

  float x, y, z; // center of mass of the node
  float mass;    // total mass of the node
  int nBodies;

  Node(Key key, bool isLeaf);
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
  Node *root;
  int leafLength;
  int nLevels; // does not count root as a level

  Octree(int maxKeyLength);
  void addChild(Node *parent, Node *child, int index, bool isLeaf, int key);
  int insert(Body body);
  void printTree(Node *node, int level);
  void buildDFT(std::vector<DFTNode> &nodes, Body *bodies);
  void traverse(Node *node, std::vector<DFTNode> &nodes);
  void setSubtreeSizes(Node *node, Body *bodies);
  void splitNode(Node *current, int index, int nLevels, int level);
};

#endif // OCTREE_H
