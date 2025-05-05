#include "octree.h"
#include <assert.h>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <string>

#define SHIFT_DIGITS 20

Node::Node(Key key, bool isLeaf) {
  this->key = key;
  this->isLeaf = isLeaf;
  this->subTreeSize = 1;
  this->whichChildren = 0;

  for (int i = 0; i < 8; ++i)
    this->children[i] = nullptr;
  this->parent = nullptr;
}

Octree::Octree(int maxKeyLength) {
  root = new Node(1, false);
  leafLength = maxKeyLength;
  nLevels = maxKeyLength / 3;
}

void Octree::addChild(Node *parent, Node *child, int index, bool isLeaf,
                      int mykey) {
  assert(parent->isLeaf == false);
  assert(mykey > 0);

  parent->children[index] = new Node(mykey, isLeaf);
  parent->children[index]->parent = parent;
  parent->whichChildren |= 1 << index;
}

void Octree::splitNode(Node *current, int index, int nLevels, int level) {
  Node *child = current->children[index];
  Key shifted = child->key >> (3 * (nLevels - level - 1));
  int childIndex = (int)shifted & 0x7;
  current->children[index] = new Node(shifted, false);
  current->children[index]->parent = current;

  child->parent = current->children[index];
  current->children[index]->children[childIndex] = child;
  current->children[index]->whichChildren |= 1 << childIndex;
}

int Octree::insert(Body body) {
  Node *current = root;
  Key key = body.key;

  int arrIndex = body.index;

  assert(key >> (leafLength - 1) == 1);  // check that all keys are prepended
  assert(nLevels * 3 == leafLength - 1); // check nlevels

  // insert internal nodes
  int level = 0;
  while (level < nLevels - 1) {
    Key shifted = key >> (3 * (nLevels - level - 1));
    int index = (int)shifted & 0x7;

    if (current->whichChildren == 0) {
      break;
    }

    if (current->children[index] == NULL) {
      addChild(current, current->children[index], index, false, shifted);

      for (int i = 0; i < 8; i++) {
        if (current->children[i] != NULL && current->children[i]->isLeaf) {
          splitNode(current, i, nLevels, level);
        }
      }
    } else {
      Node *child = current->children[index];
      if (child->isLeaf) {
        splitNode(current, index, nLevels, level);
      }
    }

    current = current->children[index];
    level++;
  }

  // insert leaf
  // assert(level == nLevels - 1);
  int index = (key >> (3 * (nLevels - level - 1))) & 0x7;

  // otherwise have to handle duplicate keys
  if (current->children[index] == NULL) {
    addChild(current, current->children[index], index, true, key);
    current->children[index]->bodies.push_back(arrIndex);

    return 1;

  } else {
    current->children[index]->bodies.push_back(arrIndex);
    if (current->children[index]->bodies.size() > 1) {
    }
    return 0;
  }
};

void Octree::printTree(Node *node, int level) {
  if (node == NULL)
    return;

  for (int i = 0; i < level; i++)
    printf("  ");

  printf("0b%s (%lu)\n", binaryString<Key>(node->key).c_str(), node->key);

  if (!node->isLeaf) {
    for (int i = 0; i < 8; i++) {
      printTree(node->children[i], level + 1);
    }
  }
};

Key getKey(Body body, int shift_digits) {
  Key key = 0;

  float x = body.x, y = body.y, z = body.z;

  // shift and trunacte
  Key shift = 1 << shift_digits;
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
  assert(pos <= (sizeof(Key) * 8 - 1));
  key += (1LL << (sizeof(Key) * 8 - 1)); // prepend 1
  printf("key: %lu\n", key);

  return key;
}

// should be a multiple of 3
int binaryLength(Key n) {
  if (n == 0)
    return 0;

  int length = (int)log2(n) + 1;
  if (length % 3 != 0)
    length += 3 - length % 3;
  assert(length % 3 == 0);
  return length;
}

void Octree::setSubtreeSizes(Node *node, Body *bodies) {
  if (node == NULL)
    return;

  float x = 0.0f, y = 0.0f, z = 0.0f;
  float mass = 0.0f;

  if (node->isLeaf) {
    node->subTreeSize = 1;

    // compute center of mass and relevant aggregate information

    if (node->bodies.size() == 0) {
      return;
    }

    for (int i = 0; i < node->bodies.size(); i++) {
      mass += bodies[node->bodies[i]].m;
      x += bodies[node->bodies[i]].x;
      y += bodies[node->bodies[i]].y;
      z += bodies[node->bodies[i]].z;
    }
    node->x = x / node->bodies.size();
    node->y = y / node->bodies.size();
    node->z = z / node->bodies.size();
    node->mass = mass;
    node->nBodies = node->bodies.size();
    return;
  }

  for (int i = 0; i < 8; i++) {
    if (node->whichChildren & (1 << i)) {
      setSubtreeSizes(node->children[i], bodies);
      node->subTreeSize += node->children[i]->subTreeSize;

      x += node->children[i]->x * node->children[i]->nBodies;
      y += node->children[i]->y * node->children[i]->nBodies;
      z += node->children[i]->z * node->children[i]->nBodies;
      mass += node->children[i]->mass;
      node->nBodies += node->children[i]->nBodies;
    }
  }

  if (node->nBodies > 0) {
    node->x = x / node->nBodies;
    node->y = y / node->nBodies;
    node->z = z / node->nBodies;
    node->mass = mass;
  }
}

void Octree::buildDFT(std::vector<DFTNode> &nodes, Body *bodies) {
  setSubtreeSizes(this->root, bodies);
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
  dftNode.isLeaf = node->isLeaf;
  dftNode.index = nodes.size();

  dftNode.x = node->x;
  dftNode.y = node->y;
  dftNode.z = node->z;
  dftNode.mass = node->mass;
  dftNode.nBodies = node->nBodies;

  dftNode.autorope =
      dftNode.index +
      node->subTreeSize; // if autorope is out of bounds, traversal is donexw

  dftNode.bodies = node->bodies;

  nodes.push_back(dftNode);

  if (!node->isLeaf) {
    for (int i = 0; i < 8; i++) {
      if (node->whichChildren & (1 << i))
        traverse(node->children[i], nodes);
    }
  }
}
