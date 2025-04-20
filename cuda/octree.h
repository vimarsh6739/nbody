#ifndef OCTREE_BUILDER_H
#define OCTREE_BUILDER_H

#include <vector>
#include <string>

using Key = long long;

struct Body {
    float x,y,z;
    float vx,vy,vz;
    float m;
    Key   key;
    int   index;
};

struct DFTNode {
    Key   key;
    bool  isLeaf;
    int   autorope;
    float x,y,z;
    float mass;
    int   nBodies;
    int   singleBodyIndex;
    float size;
};

class OctreeNode {
public:
    Key key;
    bool isLeaf;
    OctreeNode* parent;
    OctreeNode* children[8]{};
    char whichChildren = 0;
    std::vector<int> bodies;
    int   subTreeSize = 1;
    float x=0,y=0,z=0,mass=0;
    int   nBodies=0;
    int   singleBodyIndex=-1;

    OctreeNode(Key k,bool leaf,OctreeNode* p=nullptr);
    ~OctreeNode() = default;
};

class OctreeBuilder {
private:
    OctreeNode* root = nullptr;
    int keyBits;
    int nLevels;
    int totalBodies;

    void deleteRecursive(OctreeNode* node);
    OctreeNode* addChild(OctreeNode* parentNode, int index, Key childKey, bool isLeaf);
    OctreeNode* splitLeaf(OctreeNode* parentNode, int index, int parentLevel);
    void aggregate(OctreeNode* node, const Body* allBodies);
    void buildDFTRecursive(OctreeNode* node, std::vector<DFTNode>& dftArray, int depth);

public:
    OctreeBuilder(int keyLength, int numBodies);
    ~OctreeBuilder();
    void insert(Body& body);
    void buildDFT(std::vector<DFTNode>& dftOutput, const Body* allBodies);
};

class PhiloxEngine;

Key mortonNoPrepend(const Body& b);
int binaryLength(Key k);
void randomInit(Body* b, int n, PhiloxEngine& rng);

#endif
