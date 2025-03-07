#include "octree.h"
#include <cstdlib>
#include <math.h>
#include <stdio.h>

#include <assert.h>
#include <math.h>

Node::Node(int key, bool isLeaf) {
  this->key = key;
  this->isLeaf = isLeaf;
}

Octree::Octree(int maxKeyLength) {
  root = new Node(1, false);
  leafLength = maxKeyLength;
  nLevels = maxKeyLength / 3;
  printf("Octree with %d levels\n", nLevels);
}

void Octree::addChild(Node *parent, Node *child, int index, bool isLeaf,
                      int mykey) {
  assert(parent->isLeaf == false);
  assert(mykey > 0);
  parent->children[index] = new Node(mykey, false);
  parent->children[index]->parent = parent;
  parent->whichChildren |= 1 << index;
}

void Octree::insert(Body body) {
  Node *current = root;
  int key = body.key;

  assert(key >> (leafLength - 1) == 1);  // check that all keys are prepended
  assert(nLevels * 3 == leafLength - 1); // check nlevels

  // insert internal nodes
  int level = 0;
  while (level < nLevels - 1) {
    Key shifted = key >> (3 * (nLevels - level - 1));
    int index = shifted & 0x7;

    if (current->children[index] == NULL) {
      addChild(current, current->children[index], index, false, shifted);
    }

    current = current->children[index];
    level++;
  }

  // insert leaf
  assert(level == nLevels - 1);
  int index = key & 0x7;

  // otherwise have to handle duplicate keys
  assert(current->children[index] == NULL);
  addChild(current, current->children[index], index, true, key);
};

void Octree::printTree(Node *node, int level) {
  if (node == NULL)
    return;

  for (int i = 0; i < level; i++)
    printf("  ");

  printf("0b%s\n", binaryString(node->key).c_str());

  if (!node->isLeaf) {
    for (int i = 0; i < 8; i++) {
      printTree(node->children[i], level + 1);
    }
  }
};

std::string binaryString(Key k) {
  std::string s = "";
  if (k > 1)
    s = binaryString(k / 2);

  return s + std::to_string(k % 2);
}

int getKeyNoPrepend(Body body) {
  int key = 0;

  // interleave x, y, z bits
  int x = body.x;
  int y = body.y;
  int z = body.z;

  int i = 0;
  while (x || y || z) {
    key |= (body.x & (1 << i)) << i;
    key |= (body.y & (1 << i)) << (i + 1);
    key |= (body.z & (1 << i)) << (i + 2);

    x = x >> 1;
    y = y >> 1;
    z = z >> 1;
    i += 3;
  }

  // assert no key overflow
  assert(i <= sizeof(Key) * 8 - 1);

  return key;
}

// should be a multiple of 3
int binaryLength(int n) {
  if (n == 0)
    return 0;

  int length = (int)log2(n) + 1;
  if (length % 3 != 0)
    length += 3 - length % 3;
  assert(length % 3 == 0);
  return length;
}