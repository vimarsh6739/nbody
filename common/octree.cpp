#include "octree.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

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

Node::Node(Key key, bool isLeaf) {
  this->key = key;
  this->isLeaf = isLeaf;
  this->subTreeSize = 1;
  this->whichChildren = 0;

  for (int i = 0; i < 8; ++i)
    this->children[i] = nullptr;
  this->parent = nullptr;
}

Node::~Node() {
  for (int i = 0; i < 8; ++i) {
    delete children[i];
    children[i] = nullptr;
  }
}

Octree::Octree(int maxKeyLength, int resolution)
    : resolution(resolution), root(nullptr), leafLength(maxKeyLength),
      nLevels(maxKeyLength / 3) {
  root = new Node(1, false);
}

Octree::~Octree() {
  delete root;
  root = nullptr;
}

bool Octree::isLeaf(Node *node) {
  return (node->key & (1LL << (3 * resolution))) != 0;
}

void Octree::addChild(Node *parent, int index, bool isLeaf, Key mykey) {
  if (!isLeaf) {
    mykey = (parent->key << 3LL) | index;
  }

  assert(parent->isLeaf == false);
  assert(mykey > 0);

  parent->children[index] = new Node(mykey, isLeaf);
  parent->children[index]->parent = parent;
  parent->whichChildren |= 1 << index;
}

void Octree::splitNode(Node *current, int index, int nLevels, int level) {
  // pushes the child of current one level deeper into the tree
  Node *child = current->children[index]; // child needs to be split

  int currentlevel = binaryLength(current->key) / 3;
  int newlevel = currentlevel + 1;
  int childlevel = binaryLength(child->key) / 3;

  assert(newlevel < nLevels);
  assert(newlevel > currentlevel);
  assert(childlevel > newlevel);

  // compute the key for the new node
  Key newkey = child->key >> (3 * (nLevels - newlevel));

  // compute index for child as a child of the new node
  int childIndex = (child->key >> (3 * (nLevels - (newlevel + 1)))) & 0x7;

  Node *newParent = new Node(newkey, false);

  current->children[index] = newParent;
  newParent->parent = current;

  child->parent = newParent;
  newParent->children[childIndex] = child;
  newParent->whichChildren |= 1 << childIndex;
}

int Octree::getOctantIdx(Key key, int level) {
  Key shifted = key >> (3 * (nLevels - level - 1));
  return ((int)shifted & 0x7);
}

// recursive insert into octree
static void insertRec(Node *curr, Body &b) {}

int Octree::insert(Body body) {

  Node *current = root;
  Key key = body.key;

  int arrIndex = body.index;

  assert(key >> (leafLength - 1) == 1);  // check that all keys are prepended
  assert(nLevels * 3 == leafLength - 1); // check nlevels

  // insert internal nodes
  int level = 0;
  while (level < nLevels) {

    if (current->whichChildren == 0) {
      break;
    }

    int index = this->getOctantIdx(key, level);

    if (current->children[index] == nullptr) {
      // new child at this leaf must be created
      addChild(current, index, false, key);

      // if other leaves exist at this level, split them.
      // not sure if this isnecessary
      for (int i = 0; i < 8; i++) {
        if (current->children[i] != nullptr && current->children[i]->isLeaf) {
          splitNode(current, i, nLevels, level);
        }
      }
    } else {
      // child already exists. either keep traversing, or split
      Node *child = current->children[index];
      if (child->isLeaf) {
        if (child->key == key) {
          // duplicate key, chain body into bucket
          break;
        }
        splitNode(current, index, nLevels, level);
      }
    }

    current = current->children[index];
    level = binaryLength(current->key) / 3;
  }

  // insert leaf
  int index = (key >> (3 * (nLevels - level - 1))) & 0x7;
  int newLeafCreated = 0;

  // otherwise have to handle duplicate keys
  if (current->children[index] == nullptr) {
    addChild(current, index, true, key);
    assert(current->children[index]->key == key);
    current->children[index]->bodyIdx.push_back(arrIndex);
    newLeafCreated = 1;
  } else {
    assert(current->children[index]->key == key);
    current->children[index]->bodyIdx.push_back(arrIndex);
    if (current->children[index]->bodyIdx.size() > 1) {
    }
  }

  return newLeafCreated;
};

void Octree::printTree(Node *node, int level) {
  if (node == NULL)
    return;

  for (int i = 0; i < level; i++)
    printf("  ");

  if (node->isLeaf) {
    printf("0b%s (%lu) (%d)\n", binaryString<Key>(node->key).c_str(), node->key,
           node->nBodies);
  } else {
    printf("0b%s (%lu)\n", binaryString<Key>(node->key).c_str(), node->key);
  }

  if (!node->isLeaf) {
    for (int i = 0; i < 8; i++) {
      printTree(node->children[i], level + 1);
    }
  }
};

void Octree::setSubtreeSizes(Node *node, Body *bodies) {
  if (node == NULL)
    return;

  float x = 0.0f, y = 0.0f, z = 0.0f;
  float mass = 0.0f;

  node->nLeaves = 0;

  if (node->isLeaf) {
    node->subTreeSize = 1;

    node->nLeaves = node->bodyIdx.size();

    // compute center of mass and relevant aggregate information
    if (node->bodyIdx.size() == 0) {
      return;
    }

    for (int i = 0; i < node->bodyIdx.size(); i++) {
      mass += bodies[node->bodyIdx[i]].m;
      x += bodies[node->bodyIdx[i]].x;
      y += bodies[node->bodyIdx[i]].y;
      z += bodies[node->bodyIdx[i]].z;
    }
    node->x = x / node->bodyIdx.size();
    node->y = y / node->bodyIdx.size();
    node->z = z / node->bodyIdx.size();
    node->mass = mass;
    node->nBodies = node->bodyIdx.size();
    return;
  }

  for (int i = 0; i < 8; i++) {
    if (node->whichChildren & (1 << i)) {
      setSubtreeSizes(node->children[i], bodies);
      node->subTreeSize += node->children[i]->subTreeSize;
      node->nLeaves += node->children[i]->nLeaves;

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

  dftNode.bodies = node->bodyIdx;

  nodes.push_back(dftNode);

  if (!node->isLeaf) {
    for (int i = 0; i < 8; i++) {
      if (node->whichChildren & (1 << i))
        traverse(node->children[i], nodes);
    }
  }
}
