#include <cmath>
#include <gtest/gtest.h>

#include "main.h"

// Test getKeyfunction
TEST(OctreeTest, GetKey) {
  // Create a test body
  Body body = {0.5f, 0.25f, 0.125f, 0.0f, 0.0f, 0.0f, 1.0f, 0, 0};

  // Get key
  Key key = getKey(body, 4);

  // Since the function interleaves bits from x, y, z, we can test
  // specific properties of the resulting key
  ASSERT_GT(key, 0);

  // Create another test case with different values
  Body body2 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0, 0};
  Key key2 = getKey(body2, 20);
  ASSERT_EQ(key2, 1LL << 60);

  // Test with coordinates at the max (1.0): this actually should overflow our
  // resolution of 15 bits (since 1 << 15 actually has 16 bits)
  Body body3 = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0, 0};
  Key key3 = getKey(body3, 15);

  Key expectedKey3 = (7LL << (Key)(15 * 3)) + (1LL << (45));
  ASSERT_EQ(key3,
            expectedKey3); //+ 1LL << (Key)(sizeof(Key) * 8 - 1));

  // Different keys for different bodies
  ASSERT_NE(key, key2);
  ASSERT_NE(key, key3);
  ASSERT_NE(key2, key3);
}

// Test binaryLength function
TEST(OctreeTest, BinaryLength) {

  // Test case 1: Zero
  ASSERT_EQ(binaryLength(0), 0);

  // Test case 2: Powers of 2
  ASSERT_EQ(binaryLength(1), 1); // 1 in binary is 1, length is 1
  ASSERT_EQ(binaryLength(2), 2); // 2 in binary is 10, length is 2
  ASSERT_EQ(binaryLength(4), 3); // 4 in binary is 100, length is 3
  ASSERT_EQ(binaryLength(8), 4); // 8 in binary is 1000, length is 4

  // Test case 3: Random numbers
  ASSERT_EQ(binaryLength(7), 3);  // 7 in binary is 111, length is 3
  ASSERT_EQ(binaryLength(15), 4); // 15 in binary is 1111, length is 4
  ASSERT_EQ(binaryLength(16), 5); // 16 in binary is 10000, length is 5

  // Test case 4: Large numbers
  ASSERT_EQ(binaryLength(1 << 20), 21); // 1 followed by 20 zeros, length is 21
}

// Test binaryString function
TEST(OctreeTest, BinaryString) {
  ASSERT_EQ(binaryString((short)0), "0");
  ASSERT_EQ(binaryString((short)1), "1");
  ASSERT_EQ(binaryString((short)2), "10");
  ASSERT_EQ(binaryString((short)3), "11");
  ASSERT_EQ(binaryString((short)4), "100");
  ASSERT_EQ(binaryString((short)5), "101");
  ASSERT_EQ(binaryString((short)7), "111");
  ASSERT_EQ(binaryString((short)8), "1000");
  ASSERT_EQ(binaryString((short)10), "1010");
  ASSERT_EQ(binaryString((short)15), "1111");
  ASSERT_EQ(binaryString((short)16), "10000");
}

// Test Node constructor
TEST(OctreeTest, NodeConstructor) {
  Node node(42, true);

  // Check initial values
  ASSERT_EQ(node.key, 42);
  ASSERT_TRUE(node.isLeaf);
  ASSERT_EQ(node.subTreeSize, 1);
  ASSERT_EQ(node.maskChildren, 0);
  ASSERT_EQ(node.parent, nullptr);

  // Check that all children are nullptr
  for (int i = 0; i < 8; ++i) {
    ASSERT_EQ(node.children[i], nullptr);
  }
}

// Test Octree constructor
TEST(OctreeTest, OctreeConstructor) {
  Octree octree(7);

  // Check initial values
  ASSERT_EQ(octree.resolution, 7); // 21 / 3 = 7
  ASSERT_EQ(octree.root, nullptr);
  ASSERT_TRUE(octree.isEmpty(octree.root));
}

// Test Octree insertion
TEST(OctreeTest, OctreeInsert) {
  // GTEST_SKIP() << "Skipping octree insertion";
  Octree octree(16);

  // Create a test body with a properly formatted key
  // The key must have a 1 in the most significant bit (prepended)
  // which means bit position leafLength - 1 (20 in this case)
  Body body = {0.5f, 0.25f, 0.125f, 0.0f, 0.0f, 0.0f, 1.0f, 0, 0};

  // For testing, we need to set a valid key (with 1 prepended)
  Key baseKey = getKey(body, 16);
  printf("KEY: %lu, KEY LENGTH: %d\n", baseKey, binaryLength(baseKey));
  body.key = baseKey;

  // Insert the body and check the result
  int result = octree.insert(body);
  ASSERT_EQ(result, 1); // Should return 1 for successful insertion
  
  // Insert the same body again to test handling of duplicates
  int resultDup = octree.insert(body);

  ASSERT_EQ(resultDup, 0); // Should return 0 for duplicate
}

// Test DFT error
TEST(OctreeTest, DFTError) {
  GTEST_SKIP() << "Skipping end to end test";
  Octree octree(16);

  // Create a test body with a properly formatted key
  // The key must have a 1 in the most significant bit (prepended)
  // which means bit position leafLength - 1 (20 in this case)
  int nBodies = 100;
  Body *particles = new Body[nBodies];

  // Initialize positions, velocities, and masses (7 values per body)
  PhiloxEngine rng(2025);
  randomizeBodies(rng, particles, nBodies);

  // make a copy of bodies for accuracy check
  Body *particles_cp = new Body[nBodies];
  for (int i = 0; i < nBodies; i++) {
    particles_cp[i] = particles[i];
  }

  int nIters = 1;
  float dt = 10;
  nbodyIterate(particles, dt, nBodies, nIters);
  float err = checkAccuracy(particles, particles_cp, nBodies, nIters);

  ASSERT_LE(err, 1e-9); // Should return 0 for duplicate
}
