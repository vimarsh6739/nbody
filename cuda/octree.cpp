#include "octree.h"

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include "philox_engine.h"

OctreeNode::OctreeNode(Key k,bool leaf,OctreeNode* p)
    : key(k), isLeaf(leaf), parent(p) {
}

OctreeBuilder::OctreeBuilder(int keyLength, int numBodies)
    : keyBits(keyLength), totalBodies(numBodies) {
    int actualKeyBits = keyBits > 0 ? keyBits - 1 : 0;
    nLevels = actualKeyBits / 3;
    if (actualKeyBits > 0 && actualKeyBits % 3 != 0) {
    }
    if (nLevels <= 0 && totalBodies > 0) {
        nLevels = 1;
    }
    root = new OctreeNode(1, false);
}

OctreeBuilder::~OctreeBuilder() {
    deleteRecursive(root);
}

void OctreeBuilder::deleteRecursive(OctreeNode* node) {
    if(!node) return;
    if(!node->isLeaf) {
        for(int i=0; i<8; ++i) {
            if(node->whichChildren & (1<<i)) {
                deleteRecursive(node->children[i]);
            }
        }
    }
    delete node;
}

OctreeNode* OctreeBuilder::addChild(OctreeNode* parentNode, int index, Key childKey, bool isLeaf) {
    if (index < 0 || index >= 8) { throw std::out_of_range("addChild: index out of range"); }
    if (parentNode->isLeaf) { throw std::logic_error("addChild: cannot add child to a leaf"); }
    if (parentNode->children[index] != nullptr) {
         std::cerr << "Warning: addChild overwriting child at index " << index << std::endl;
         deleteRecursive(parentNode->children[index]);
    }
    auto* child = new OctreeNode(childKey, isLeaf, parentNode);
    parentNode->children[index] = child;
    parentNode->whichChildren |= (1 << index);
    return child;
}

OctreeNode* OctreeBuilder::splitLeaf(OctreeNode* parentNode, int index, int parentLevel) {
    if (index < 0 || index >= 8) { throw std::out_of_range("splitLeaf: index out of range"); }
    OctreeNode* originalLeaf = parentNode->children[index];
    if (!originalLeaf || !originalLeaf->isLeaf) { throw std::logic_error("splitLeaf: target is not a leaf"); }

    Key newNodeKey = originalLeaf->key;
    OctreeNode* newNode = new OctreeNode(newNodeKey, false, parentNode);

    parentNode->children[index] = newNode;

    int childLevel = parentLevel + 1;
    int shift = keyBits - 3 * (childLevel + 1);
    if (shift < 0) shift = 0;
    int childIndexForOldLeaf = (originalLeaf->key >> shift) & 7;

     if (childIndexForOldLeaf < 0 || childIndexForOldLeaf >= 8) {
        std::cerr << "Error: splitLeaf invalid child index " << childIndexForOldLeaf << std::endl;
         childIndexForOldLeaf = 0;
     }

    newNode->children[childIndexForOldLeaf] = originalLeaf;
    newNode->whichChildren |= (1 << childIndexForOldLeaf);
    originalLeaf->parent = newNode;

    return newNode;
}

void OctreeBuilder::aggregate(OctreeNode* node, const Body* allBodies) {
    if (!node) return;
    node->subTreeSize = 1; node->mass = 0; node->x = 0; node->y = 0; node->z = 0; node->nBodies = 0; node->singleBodyIndex = -1;
    if(node->isLeaf) {
        float com_x = 0, com_y = 0, com_z = 0;
        for(int idx : node->bodies) { if (idx < 0 || idx >= totalBodies) { continue; } const auto& body = allBodies[idx]; node->mass += body.m; com_x += body.x * body.m; com_y += body.y * body.m; com_z += body.z * body.m; }
        node->nBodies = node->bodies.size();
        if(node->mass > 1e-20f){ node->x=com_x/node->mass; node->y=com_y/node->mass; node->z=com_z/node->mass; }
        else if (node->nBodies > 0) { node->x=allBodies[node->bodies[0]].x; node->y=allBodies[node->bodies[0]].y; node->z=allBodies[node->bodies[0]].z; }
        if(node->nBodies == 1) { node->singleBodyIndex = node->bodies[0]; }
    } else {
        float total_children_mass = 0; float com_x_num = 0, com_y_num = 0, com_z_num = 0;
        for(int i=0; i<8; ++i) {
            if (node->whichChildren & (1 << i)) { // Check if child exists
                OctreeNode* child = node->children[i];
                if (!child) { // Safety check
                     std::cerr << "Warning: aggregate null child ptr for index " << i << " (Node key: " << node->key << ").\n";
                     continue;
                }
                aggregate(child, allBodies);
                node->subTreeSize += child->subTreeSize;
                node->nBodies += child->nBodies;
                if (child->mass > -1e-20f) {
                    total_children_mass += child->mass;
                    com_x_num += child->x * child->mass;
                    com_y_num += child->y * child->mass;
                    com_z_num += child->z * child->mass;
                }
            }
        }
        node->mass = total_children_mass;
        if(node->mass > 1e-20f){ node->x=com_x_num/node->mass; node->y=com_y_num/node->mass; node->z=com_z_num/node->mass; }
    }
}

void OctreeBuilder::buildDFTRecursive(OctreeNode* node, std::vector<DFTNode>& dftArray, int depth) {
     if (!node) return;
    DFTNode dftNode{}; dftNode.key = node->key; dftNode.isLeaf = node->isLeaf; dftNode.autorope = dftArray.size() + node->subTreeSize;
    dftNode.x = node->x; dftNode.y = node->y; dftNode.z = node->z; dftNode.mass = node->mass; dftNode.nBodies = node->nBodies; dftNode.singleBodyIndex = node->singleBodyIndex;
    dftNode.size = powf(0.5f, depth); dftArray.push_back(dftNode);
    if (!node->isLeaf) { for (int i = 0; i < 8; ++i) { if (node->whichChildren & (1 << i)) { if (node->children[i]) { buildDFTRecursive(node->children[i], dftArray, depth + 1); } } } }
}

void OctreeBuilder::insert(Body& body) {
    OctreeNode* current = root; Key key = body.key;
    if (!root) throw std::logic_error("insert: Root is null");
    if (nLevels <= 0) { root->bodies.push_back(body.index); return; }

    for (int level = 0; level < nLevels; ++level) {
        int shift = keyBits - 3 * (level + 1); if (shift < 0) shift = 0;
        int index = (key >> shift) & 7;
        OctreeNode* childNode = current->children[index];

        if (childNode == nullptr) {
            childNode = addChild(current, index, key, (level == nLevels - 1));
        } else if (childNode->isLeaf && level < nLevels - 1) {
             childNode = splitLeaf(current, index, level);
        }

        current = childNode;

        if (current->isLeaf) {
            current->bodies.push_back(body.index);
            return;
        }
    }
     std::cerr << "Error: insert loop finished for body " << body.index << " (key " << key << ")." << std::endl;
}

void OctreeBuilder::buildDFT(std::vector<DFTNode>& dftOutput, const Body* allBodies) {
    if (!root) return;
    aggregate(root, allBodies);
    dftOutput.clear();
    if (root->subTreeSize > 0) {
        dftOutput.reserve(root->subTreeSize);
    }
    buildDFTRecursive(root, dftOutput, 0);
}


Key mortonNoPrepend(const Body& b) {
    const int shift_digits = 20;
    float x=std::max(0.0f,std::min(b.x,0.9999999f));
    float y=std::max(0.0f,std::min(b.y,0.9999999f));
    float z=std::max(0.0f,std::min(b.z,0.9999999f));
    Key scale = (Key)1 << shift_digits;
    Key xi=(Key)(x*scale), yi=(Key)(y*scale), zi=(Key)(z*scale);
    Key key=0;
    for(int i=0; i<shift_digits; ++i) {
        key |= (((xi >> i) & 1) << (3*i + 0));
        key |= (((yi >> i) & 1) << (3*i + 1));
        key |= (((zi >> i) & 1) << (3*i + 2));
    }
    return key;
}

int binaryLength(Key k) {
    if (k <= 0) return 0;
    int len = 0; Key temp = k;
    while (temp > 0) { temp >>= 1; ++len; }
    if (len == 0) return 0;
    if (len % 3 != 0) { len += 3 - (len % 3); }
    return len;
}

void randomInit(Body* b, int n, PhiloxEngine& rng) {
     std::uniform_real_distribution<float> U01(0.0f, 1.0f);
     std::uniform_real_distribution<float> V(-0.005f, 0.005f);
     std::uniform_real_distribution<float> M(0.1f, 1.0f);
     for(int i=0; i<n; ++i) {
        b[i].x=std::max(1e-6f, U01(rng));
        b[i].y=std::max(1e-6f, U01(rng));
        b[i].z=std::max(1e-6f, U01(rng));
        b[i].vx=V(rng); b[i].vy=V(rng); b[i].vz=V(rng);
        b[i].m=M(rng); b[i].index=i;
     }
}
