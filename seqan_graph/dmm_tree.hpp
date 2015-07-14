#ifndef __DMM_TREE_HPP__
#define __DMM_TREE_HPP__

/**
 * This class implements a dynamic min-max tree as it is described by
 * Navarro and Sadakane here:
 * http://dl.acm.org/citation.cfm?id=2601073
 *
 * We take some shortcuts for now for testing feasibility that can be improved
 * later on.
 */

#include <assert.h>

#include <vector>
#include <list>

using namespace std;

// forward declarations
class TreeNode;

/**
 * This class encapsulates a single storage cell that 
 * can contain parts of up to three different storage segments
 * each of length between L and 2L.
 */
class ListNode {

    private:
        size_t L;
        bool *cell;
        TreeNode *node, *node2;

    public:
        ListNode(size_t _L) : L(_L) {
            cell = new bool[2*L]; 
        }

        bool get(size_t idx) {
            return cell[idx];
        }

        void put(bool val, size_t idx) {
            cell[idx] = val;
        }

        void setNode(size_t offset, TreeNode* _node) {
            if (offset < L)
                node = _node;
            else
                node2 = _node;
        }

        TreeNode* getNode(size_t offset) {
            return (offset < L) ? node : node2;
        }

};

class TreeNode {

    private:
        size_t i;
        size_t j;
        size_t segLen;
        size_t offset;
        TreeNode *parent, *right, *left;
        list<ListNode*>::iterator cell;
        
        bool is_red;

    public:
        TreeNode() : i(0), j(0), parent(NULL), right(NULL), left(NULL), is_red(false) {}

        TreeNode(size_t _i, size_t _j, TreeNode* _parent) : 
            i(_i), 
            j(_j), 
            parent(_parent), 
            right(NULL), 
            left(NULL), 
            is_red(true) {}

        TreeNode(size_t _i, size_t _j, TreeNode* _parent, TreeNode* _left, TreeNode* _right) : 
            i(_i), 
            j(_j), 
            parent(_parent), 
            right(_right), 
            left(_left), 
            is_red(true) {}

        /**
         * Here start the setters ...
         */
        void setLeft(TreeNode *node) {
            left = node;
        }

        void setRight(TreeNode *node) {
            right = node;
        }

        void setParent(TreeNode *node) {
            parent = node;
        }

        void setOffset(size_t _offset) {
            offset = _offset;
        }

        void setSegLen(size_t _segLen) {
            segLen = _segLen;
        }

        void setCell(list<ListNode*>::iterator _cell) {
            cell = _cell;
        }

        void setRed() {
            is_red = true;
        }

        void setBlack() {
            is_red = false;
        }

        /**
         * Here start the getters ...
         */
        TreeNode* getLeft() {
            return left;
        }

        TreeNode* getRight() {
            return right;
        }

        TreeNode* getParent() {
            return parent;
        }

        TreeNode* getGrandparent() {
            return (parent == NULL) ? NULL : parent->getParent();
        }

        TreeNode* getUncle() {
            TreeNode* node = this->getGrandparent();
            return (node == NULL) ? NULL : (parent == node->getLeft()) ? node->getRight() : node->getLeft();
        }

        size_t getOffset() {
            return offset;
        }

        size_t getSegLen() {
            return segLen;
        }

        bool isRed() {
            return is_red;
        }

        size_t getUpper() {
            return j;
        }

        size_t getLower() {
            return i;
        }

        list<ListNode*>::iterator getCell() {
            return cell;
        }

};

class SegmentStore {

    private:
        vector<list<ListNode*> > listVector;
        vector<size_t> header;
        size_t L;

    public:
        SegmentStore(size_t _L) : L(_L) {
            // create list of segments
            for (size_t i = 0; i <= L; i++) {
                list<ListNode*> tmp;
                tmp.push_front(new ListNode(L));
                listVector.push_back(tmp);
                header.push_back(2*L);
            }
        }

        /** 
         * Insert a full new segment that is taken from P[s, t). 
         */
        void insertSegment(list<bool> P, size_t s, size_t t, TreeNode *node) {

            // get segment length --> list we want to add to (P has size between L and 2L)
            size_t segLen = (t - s); 
            size_t segLenIdx = segLen - L;

            // header offset and cell iterator
            size_t sidx = header.at(segLenIdx);
            list<ListNode*>::iterator it = listVector.at(segLenIdx).begin();

            // check if there is enough room left in the current cell or wether 
            // we need to open a new one
            if (sidx < segLen) {
                listVector.at(segLenIdx).push_front(new ListNode(L));
                sidx = 2*L - (sidx - segLen);
                it = listVector.at(segLenIdx).begin();
                header.at(segLenIdx) = sidx;
            }

            // tell the node where it can find its segment
            node->setOffset(sidx);
            node->setCell(it);
            node->setSegLen(t - s);
            // tell the segment about its node
            listVector.at(segLenIdx).front()->setNode(sidx, node);
            // iterate over P and write values
            list<bool>::iterator itP = P.begin();
            advance(it, s);
            for (size_t i = s; i < t; i++, itP++) {
                (*it)->put(*itP, sidx++);
                if (sidx >= 2*L) {
                    it++;
                    sidx = 0;
                }
            }
        }

        /** 
         * Insert a value at idx into a segment of length segLen in cell cellIdx thereby
         * moving the segment around.
         */
        void insertValue(bool value, size_t idx, TreeNode *node) {

            /**
             * At first we move the segment from list i to list i+1 thereby
             * integrating the new value at position idx.
             */

            // get the info the node has about its segment
            size_t idxSource = node->getOffset();
            size_t segLen = node->getSegLen();
            size_t segLenIdx = segLen - L;
            list<ListNode*>::iterator itSource = node->getCell();

            // save idx and offset for later
            size_t idxOrig = idxSource;
            list<ListNode*>::iterator itOrig = itSource;

            // get target location
            list<ListNode*>::iterator itTarget = listVector.at(segLenIdx + 1).begin();
            size_t idxTarget = header.at(segLenIdx + 1);

            // check if there is enough space for the segment
            if (idxTarget < segLen + 1) {
                listVector.at(segLenIdx + 1).push_front(new ListNode(L));
                idxTarget = 2*L - (idxTarget - segLen - 1); 
                itTarget = listVector.at(segLenIdx + 1).begin();
            } else {
                idxTarget -= segLen + 1;
            }
            header.at(segLenIdx + 1) = idxTarget;

            // tell the node the new position of its segment
            node->setOffset(idxTarget);
            node->setCell(itTarget);
            node->setSegLen(segLen + 1);

            // tell the segment about its new node
            listVector.at(segLenIdx + 1).front()->setNode(idxTarget, node);
            
            // move data between locations
            size_t moved = 0;
            while (moved < segLen) {
                move(*itSource, idxSource++, *itTarget, idxTarget++);
                moved++;
                if (idxSource == 2*L) {
                    itSource++;
                    idxSource = 0;
                }
                if (idxTarget == 2*L) {
                    itTarget++;
                    idxTarget = 0;
                }
                if (moved == idx)
                    (*itTarget)->put(value, idxTarget++);
                if (idxTarget == 2*L) {
                    itTarget++;
                    idxTarget = 0;
                }
            }
            

            /**
             * As a second step, we move the previous header segment of list i
             * to the position that the segment we just moved has freed up.
             * The node is now the node that resided at the head of the list!
             */
            moveHeaderSegment(segLen, segLenIdx, itOrig, idxOrig);

        }

        /**
         * This function take the header segment of list idx and
         * moves it to cell itTarget with offset idxTarget.
         */
        void moveHeaderSegment(size_t segLen, size_t segLenIdx, list<ListNode*>::iterator itTarget, size_t idxTarget) {

            bool freeCell = false;
            size_t moved = 0;

            // get source information of header segment
            TreeNode *node = listVector.at(segLenIdx).front()->getNode(header.at(segLenIdx));
            size_t idxSource = node->getOffset();
            list<ListNode*>::iterator itSource = listVector.at(segLenIdx).begin();

            // update node information
            node->setOffset(idxTarget);
            node->setCell(itTarget);
            (*itTarget)->setNode(idxTarget, node);

            // we do not check if there is enough room for the data here!
            // we assume this is done by the caller.
            assert(segLen == node->getSegLen());

            // move data
            while (moved < segLen) {
                move(*itSource, idxSource++, *itTarget, idxTarget++);
                moved++;
                if (idxSource == 2*L) {
                    itSource++;
                    idxSource = 0;
                    freeCell = true;
                }
                if (idxTarget == 2*L) {
                    itTarget++;
                    idxTarget = 0;
                }
            }
            // update header info
            header.at(segLenIdx) = idxSource;

            // clean up freed cell in front
            if (freeCell) {
                delete listVector.at(segLenIdx).front();
                listVector.at(segLenIdx).pop_front();
            }

        }


        /**
         * Here we take the segment of node and split it into two new segments.
         * These two segments will then be associated to two new nodes node1 and node2.
         */
        void splitSegment(TreeNode* node, TreeNode* node1, TreeNode* node2) {
            
            // extract information from the calling node
            size_t idxSource = node->getOffset();
            size_t segLen = node->getSegLen();
            size_t segLenIdx = segLen - L;
            list<ListNode*>::iterator itSource = node->getCell();

            // save information of calling nade to later replace it
            size_t idxOrig = idxSource;
            list<ListNode*>::iterator itOrig = itSource;

            list<ListNode*>::iterator itTarget = listVector.at(0).begin();
            size_t idxTarget = header.at(0);

            // make sure this is the correct call
            assert(segLen == 2*L);
            assert(segLenIdx == L);

            // check if there is enough space for the two segments
            if (idxTarget < segLen) {
                listVector.at(0).push_front(new ListNode(L));
                idxTarget = 2*L - (idxTarget - segLen);
                itTarget = listVector.at(0).begin();
            } else {
                idxTarget -= segLen;
            }
            header.at(0) = idxTarget;

            // add information to first new node
            node1->setOffset(idxTarget);
            node1->setSegLen(segLen / 2);
            node1->setCell(itTarget);
            (*itTarget)->setNode(idxTarget, node1);

            // move first half of data
            size_t moved = 0;
            while (moved < segLen / 2) {
                move(*itSource, idxSource++, *itTarget, idxTarget++);
                if (idxSource == 2*L) {
                    itSource++;
                    idxSource = 0;
                }
                if (idxTarget == 2*L) {
                    itTarget++;
                    idxTarget = 0;
                }
            }

            // add information to second new node
            node2->setOffset(idxTarget);
            node2->setSegLen(segLen / 2);
            node2->setCell(itTarget);
            (*itTarget)->setNode(idxTarget, node2);

            // move second half of data
            while (moved < segLen) {
                move(*itSource, idxSource++, *itTarget, idxTarget++);
                if (idxSource == 2*L) {
                    itSource++;
                    idxSource = 0;
                }
                if (idxTarget == 2*L) {
                    itTarget++;
                    idxTarget = 0;
                }
            }

            // remove old segment that we just copied into the two new segments
            moveHeaderSegment(segLen, segLenIdx, itOrig, idxOrig);
        }


        /**
         * This helper function takes a bit at position s_idx in source and copies it to
         * position t_idx in target.
         */
        void move(ListNode *source, size_t s_idx, ListNode *target, size_t t_idx) {
            target->put(source->get(s_idx), t_idx);
        }

};



class DmmSubTree {

    private:

        // fix the min segment length to L (max segment length is 2 * L)
        size_t L;
        TreeNode *root;
        SegmentStore *seg_store;
        list<bool> P;
        size_t currentN;

    public:

        DmmSubTree(size_t _L) : L(_L) {
            // init array of linked lists storing segments
            seg_store = new SegmentStore(L);
            root = NULL;
            currentN = 0;
        }

        /**
         * Insert a new item into the tree. We handle a tree with less than
         * L values as a special case as we can not reflect that few items in 
         * our segment storage structure. Once P contains L values, we open
         * the root node an start the first segment from P.
         */
        void insert(size_t idx, bool value) {
            currentN++; 
            // if we currently have less values than our min segment size we use the array P directly
            if (currentN <= this->L) {
                if (currentN == 0)
                    P.push_front(value);
                else {
                    list<bool>::iterator it = P.begin();
                    advance(it, idx);
                    P.insert(it, value);
                }
            // we create the tree structure once we have reached minimal segment length
            } else if (currentN == this->L) {
                // create root node and insert the full segment
                root = new TreeNode(1, 4, NULL);
                root->setBlack();
                seg_store->insertSegment(P, 0, L, root);
            // start inserting at existing root node
            } else
                insert(this->root, idx, value);
        }

        size_t getLower() {
            return (root == NULL) ? 0 : root->getLower();
        }

        size_t getUpper() {
            return (root == NULL) ? 0 : root->getUpper();
        }


    private:

        /**
         * This function describes the internal red black tree logic that 
         * adds values to existing segments, generates new nodes once they become
         * necessary and balances the tree at the end of each insertion.
         */
        void insert(TreeNode *node, size_t idx, bool value) {
            // we reached a leaf
            // try to extend segment stored in the node
            if (node->getLeft() == NULL) {
                assert(node->getRight() == NULL); // no nodes with only one leaf by construction?

                // current segment has not reached max size yet, we can keep adding to it
                if (node->getSegLen() < 2*L) {
                    seg_store->insertValue(value, idx, node);

                // we need to split the node and add two children
                } else {
                    // create two children
                    TreeNode *node1 = new TreeNode(node->getLower(), node->getUpper() / 2, node);
                    TreeNode *node2 = new TreeNode(node->getUpper() / 2 + 1, node->getUpper(), node);
                    node->setCell(NULL);

                    // add the first child 
                    node->setLeft(node1);
                    rebalance(node1);

                    // the position of the second child depends on any rotations
                    // initiated by the previous child
                    if (node->getRight() != NULL) {
                        TreeNode *parent = node->getRight();
                        assert(parent->getLeft() == NULL);
                        parent->setLeft(node2);
                        node2->setParent(parent);
                    } else {
                        node->setRight(node2);
                        node2->setParent(node);
                    }
                    rebalance(node2);
                    
                    // initiate a segment split in the storage and link to the two new nodes
                    seg_store->splitSegment(node, node1, node2);
                }

                // finally update the node values back to the root
                update_intervals(node1);
                update_intervals(node2);

            // decide into which child we should proceed
            } else {
                if (idx <= node->getLeft()->getUpper()) {
                    insert(node->getLeft(), idx, value);
                } else {
                    insert(node->getRight(), idx - node->getLeft()->getUpper(), value);
                }
            }
        }

        /** 
         * This function walks through the cases of a red black tree starting from the
         * node we just inserted, eventually re-balancing the tree if necessary.
         */
        void rebalance(TreeNode *node) {
            
            TreeNode *parent = node->getParent();
            // we inserted the root
            if (parent == NULL) {
                node->setBlack();    
            // tree is still valid
            } else if (!parent->isRed()) {
                return;
            } else {
                TreeNode *grandparent = node->getGrandparent();
                TreeNode *uncle = node->getUncle();
                // both parent and uncle are red
                if (uncle != NULL && uncle->isRed()) {
                    parent->setBlack();
                    uncle->setBlack();
                    grandparent->setRed();
                    rebalance(grandparent);
                // our parent is red but the uncle is black
                } else {
                    if (node == parent->getRight() && parent == grandparent->getLeft()) {
                        rotateLeft(parent);
                        // procede with the left child
                        node = node->getLeft();
                    } else if (node == parent->getLeft() && parent == grandparent->getRight()) {
                        rotateRight(parent);
                        // procede with the right child
                        node = node->getRight();
                    }
                    // we might have changed node and need to update the (grand)parent
                    grandparent = node->getGrandparent();
                    parent = node->getParent();
                    // re-color parents
                    parent->setBlack();
                    grandparent->setRed();
                    if (node == parent->getLeft())
                        rotateRight(grandparent);
                    else
                        rotateLeft(grandparent);
                }
            }

        }

        /*
         * This function will rotate round the given node. That is, the current right child 
         * will become the father of the node.
         */
        void rotateLeft(TreeNode* node) {
            // get right child
            TreeNode *rightChild = node->getRight();
            TreeNode *parent = node->getParent();
            node->setRight(rightChild->getLeft());
            // relink father to the right child
            parent->setLeft(rightChild);
            // update parent relationships
            rightChild->setParent(parent);
            // attach the node as left child to its former right child
            rightChild->setLeft(node);
            node->setParent(rightChild);
        }

        /*
         * This function will rotate around the given node and makes the node's father its
         * right child.
         */
        void rotateRight(TreeNode* node) {
            // get left child
            TreeNode *rightChild = node->getLeft();
            // replace left child with node's parent
            TreeNode *parent = node->getParent();
            // make the node's parent its right child
            node->setRight(parent);
            // update node parenthoods
            node->setParent(parent->getParent());
            parent->setParent(node);
            // reattach the node's right child to the former parent's left
            parent->setLeft(rightChild);
            
            swapUpper(node, parent);
        }

        void swapUpper(TreeNode *node1, TreeNode *node2) {
            size_t tmp = node1->getUpper();
            node1->setUpper(node2->getUpper());
            node2->setUpper(tmp);
        }

};


/**
 * This is the tree master class that will hold the three dmmTrees that
 * are needed for memory balancing when dynamically increasing / decreasing
 * the number of elements in the array and we need to adapt the chunk size L.
 */
class DmmTree {

    private:
        DmmSubTree *current, *previous, *next;

    public:
        
        // we begin with an empty tree and use P until length 4 is reached
        DmmTree() : 
            current(NULL),
            previous(NULL),
            next(NULL) {}

        // insert value into the tree
        void insert(bool value, size_t i) {
            // decide in which sub tree to insert
            if (current == NULL) {
                current = new DmmSubTree(1);
                current->insert(value, i);
            } else {
                // insert to previous
                if (i < current->getLower()) {
                    assert(previous != NULL);
                    moveUp(previous, current);
                    previous->insert(value, i); 
                    moveUp(previous, current);
                    moveUp(current, next);
                    moveUp(current, next);
                // insert into next
                } else if (i > current->getUpper()) {
                    assert(next != NULL);
                    next->insert(value, i); 
                    moveUp(previous, current);
                    moveUp(current, next);
                // insert into current
                } else {
                    moveUp(current, next);
                    moveUp(current, next);
                    current->insert(value, i);
                    moveUp(previous, current); 
                }

                // check if the "current" sub tree is empty and we need to
                // make either "previous" or "next" the new "current", create a new
                // sub tree and delete the one no longer used.
                // TODO
            }
        }

    private:

        /** This function moves the last bit of source to the 
         * beginning of target.
         */
        void moveUp(DmmSubTree *source, DmmSubTree *target) {
            // TODO
        }

        /** This function moves the first bit of source to the 
         * end of target.
         */
        void moveDown(DmmSubTree *source, DmmSubTree *target) {
            // TODO
        }

        /**
         * Quick function to get the ceiling of log2 of 
         * an unsigned int.
         */
        inline size_t ceilLog2(size_t val) {
            size_t ret = 0;
            while (val >>= 1)
                ret++;
            return ret;
        }
};
#endif
