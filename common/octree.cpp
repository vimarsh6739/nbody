#include "octree.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

/* Morton Representation */

Key getKey(Body body, int resolution) {
  Key key = 0;

  float x = body.x, y = body.y, z = body.z;

  // shift and trunacte
  Key shift = 1LL << resolution;
  Key xint = static_cast<Key>(x * shift);
  Key yint = static_cast<Key>(y * shift);
  Key zint = static_cast<Key>(z * shift);

  // Check for overflow
  if (xint < 0 || yint < 0 || zint < 0) {
    fprintf(stderr, "Error: Integer overflow detected in getKeyNoPrepend\n");
    exit(EXIT_FAILURE);
  }

  uint64_t pos = 0;
  while (xint || yint || zint) {
    // construct using LSB
    key |= (xint & 1) << pos;
    key |= (yint & 1) << (pos + 1LL);
    key |= (zint & 1) << (pos + 2LL);

    xint = xint >> 1LL;
    yint = yint >> 1LL;
    zint = zint >> 1LL;

    pos += 3;
  }

  // assert no key overflow
  assert((pos + 1) <= sizeof(Key) * 8);

  // actual key size = `(3*shift_digits+1)` bits
  // prepend a 1 to the key (at (1<<3*shift_digits)) because its a leaf
  key += (1LL << (3 * resolution));

  return key;
}

// should be a multiple of 3
int binaryLength(Key n) {
  if (n == 0)
    return 0;

  int length = (int)log2(n) + 1;
  return length;
}

/* Node functions */

Node::Node(Key key) {
  this->key = key;
  this->subTreeSize = 1;
  this->maskChildren = 0;

  for (int i = 0; i < 8; ++i)
    this->children[i] = nullptr;
}

Node::~Node() {
  for (int i = 0; i < 8; ++i) {
    delete children[i];
    children[i] = nullptr;
  }
}

/* Octree functions */

Octree::Octree(int resolution) : resolution(resolution), root(nullptr) {}

Octree::~Octree() {
  delete root;
  root = nullptr;
}

bool Octree::isEmpty(Node *&node) { return node == nullptr; }

bool Octree::isLeaf(Node *&node) {
  return (node->key & (1LL << (3 * resolution))) != 0;
}

int Octree::getOctantIndex(Key key, int level) {
  // [LEAF_BIT] [level0 : 3] [level1 : 3] [level2 : 3] ... [level{r-1} : 3]
  Key shifted = key >> (3 * (resolution - 1 - level));
  return ((int)shifted & 0x7);
}

void Octree::updateAggregateStats(Node *&node, Body &body) {
  // update centroid and total mass
  node->cx += body.x * body.m;
  node->cy += body.y * body.m;
  node->cz += body.z * body.m;
  node->tm += body.m;

  // update subTreeSize, handle duplicates
  if (!isLeaf(node) || (node->bodyIdx.size() > 0))
    node->subTreeSize++;
}

int Octree::subdivide(Node *&node, int level) {
  int subidx = this->getOctantIndex(node->key, level + 1);
  node->children[subidx] = new Node(node->key); // push down leaf
  node->maskChildren |= (1 << subidx);                // update mask
  node->key -= (1LL << (3 * resolution));             // mark as internal node
  return subidx;
}

int Octree::insert(Body &body) {

  if (isEmpty(root)) {
    root = new Node(body.key);
    return 1;
  }

  int addedLeaf = 0;
  Node *ptr = root;
  for (int level = 0; level < resolution; ++level) {
    if (isEmpty(ptr)) {
      // alloc and mark as leaf
      ptr = new Node(body.key);
      updateAggregateStats(ptr, body);
      addedLeaf = 1;
      break;
    } else if (isLeaf(ptr)) {
      if (level < resolution - 1) {
        int next = subdivide(ptr, level);
        updateAggregateStats(ptr, body);
        ptr = ptr->children[next];
      } else {
        // duplicate key
        ptr->bodyIdx.push_back(body.index);
        updateAggregateStats(ptr, body);
      }
    } else {
      updateAggregateStats(ptr, body);
      int next = this->getOctantIndex(body.key, level + 1);
      ptr = ptr->children[next];
    }
  }

  return addedLeaf;
}

void Octree::printTree(Node *&node, int level) {

  if (isEmpty(node))
    return;

  printf("%d:", level);

  for (int i = 0; i < level; i++)
    printf("  ");

  printf("0b%s (%lu) (%lu)\n", binaryString<Key>(node->key).c_str(), node->key,
         node->subTreeSize);

  if (!isLeaf(node)) {
    for (int i = 0; i < 8; i++) {
      printTree(node->children[i], level + 1);
    }
  }
};

void Octree::buildDFT(std::vector<DFTNode> &nodes, Body *bodies) {
  traverse(this->root, nodes);

  // printf("DFT: ");
  // for (int i = 0; i < nodes.size(); i++) {
  //   int autoropeNext = nodes[i].autorope;
  //   if (autoropeNext >= nodes.size())
  //     printf("%d (-1), ", nodes[i].key);
  //   else
  //     printf("%d (%d), ", nodes[i].key, nodes[autoropeNext].key);
  // }
  // printf("\n");
}

void Octree::traverse(Node *node, std::vector<DFTNode> &nodes) {
  if (node == NULL)
    return;

  DFTNode dftNode;
  dftNode.key = node->key;
  dftNode.isLeaf = isLeaf(node);
  dftNode.index = nodes.size();

  dftNode.x = node->cx;
  dftNode.y = node->cy;
  dftNode.z = node->cz;
  dftNode.mass = node->tm;

  dftNode.autorope =
      dftNode.index +
      node->subTreeSize; // if autorope is out of bounds, traversal is donexw

  dftNode.bodies = node->bodyIdx;

  nodes.push_back(dftNode);

  if (!isLeaf(node)) {
    for (int i = 0; i < 8; i++) {
      if (node->maskChildren & (1 << i))
        traverse(node->children[i], nodes);
    }
  }
}
