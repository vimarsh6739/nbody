#ifndef OCTREE_H
#define OCTREE_H

#include <string>
#include <vector>

typedef int Key;

typedef struct {
  int x, y, z;         // position must be ints for key building
  float vx, vy, vz, m; // each body now has a mass 'm'
  Key key;
} Body;

int getKeyNoPrepend(Body b);
int binaryLength(int n);
std::string binaryString(Key key);

class Node {
public:
  int key;
  bool isLeaf;

  Node *parent;
  Node *children[8];
  char whichChildren;

  int subTreeSize; // used to build autoropes

  Node(int key, bool isLeaf);
};

typedef struct DFTNode {
  int index;
  int key;
  bool isLeaf;
  int autorope;
};
class Octree {

public:
  Node *root;
  int leafLength;
  int nLevels; // does not count root as a level

  Octree(int maxKeyLength);
  void addChild(Node *parent, Node *child, int index, bool isLeaf, int key);
  void insert(Body body);
  void printTree(Node *node, int level);
  void buildDFT(std::vector<DFTNode> &nodes);
  void traverse(Node *node, std::vector<DFTNode> &nodes);
  void setSubtreeSizes(Node *node);
};

#endif // OCTREE_H
