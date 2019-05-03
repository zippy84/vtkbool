/*
Copyright 2012-2019 Ronald Römer

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef __AABB_h
#define __AABB_h

class BB {
public:
    BB (double minX, double maxX, double minY, double maxY) : minX(minX), maxX(maxX), minY(minY), maxY(maxY) {}
    BB () {}

    double GetArea () {
        return (maxX-minX)*(maxY-minY);
    }

    std::shared_ptr<BB> Merge (const BB &other) {
        return std::make_shared<BB>(std::min(minX, other.minX),
            std::max(maxX, other.maxX),
            std::min(minY, other.minY),
            std::max(maxY, other.maxY));
    }

    bool Intersects (const BB &other) const {
        return other.maxX >= minX && other.minX <= maxX
            && other.maxY >= minY && other.minY <= maxY;
    }

    double minX, maxX, minY, maxY;
};

class Obj {
public:
    virtual BB GetBB () const = 0;
};

class Line : public Obj {
public:
    Line (const Point &a, const Point &b, int grp = NO_USE) : bb{std::min(a.x, b.x),
        std::max(a.x, b.x),
        std::min(a.y, b.y),
        std::max(a.y, b.y)},
        grp(grp), pA(a), pB(b) { }

    virtual BB GetBB () const {
        return bb;
    }

    BB bb;
    int grp;

    Point pA, pB;

    friend std::ostream& operator<< (std::ostream &out, const Line &l) {
        out << "{" << l.pA << "}, {" << l.pB << "}";
        return out;
    }
};

class Node {
public:
    BB bb;
    Node (std::shared_ptr<Obj> obj) : obj(obj), left(NO_USE), right(NO_USE), parent(NO_USE) {
        if (obj) {
            bb = obj->GetBB();
        }
    }

    int left, right, parent;
    std::shared_ptr<Obj> obj;
};

class AABB {
public:
    AABB () : rootId(NO_USE) {}

    void InsertObj (std::shared_ptr<Obj> obj) {
        int nodeId = CreateNode(obj);

        // std::cout << nodes.at(nodeId).obj.use_count() << std::endl;

        if (rootId == NO_USE) {
            rootId = nodeId;
        } else {
            int _id = rootId;

            while (nodes[_id].left != NO_USE) {
                Node &curr = nodes[_id];

                BB bb = *curr.bb.Merge(nodes[nodeId].bb);

                double diff = bb.GetArea()-curr.bb.GetArea();

                BB bbA = *nodes[curr.left].bb.Merge(nodes[nodeId].bb),
                    bbB = *nodes[curr.right].bb.Merge(nodes[nodeId].bb);

                double cLeft, cRight;

                if (nodes[curr.left].left == NO_USE) {
                    cLeft = bbA.GetArea()+diff;
                } else {
                    cLeft = bbA.GetArea()-nodes[curr.left].bb.GetArea()+diff;
                }

                if (nodes[curr.right].left == NO_USE) {
                    cRight = bbB.GetArea()+diff;
                } else {
                    cRight = bbB.GetArea()-nodes[curr.right].bb.GetArea()+diff;
                }

                if (bb.GetArea() < cLeft && bb.GetArea() < cRight) {
                    break;
                }

                if (cLeft < cRight) {
                    _id = curr.left;
                } else {
                    _id = curr.right;
                }

            }

            int parA = nodes[_id].parent;
            int parB = CreateNode();

            nodes[parB].bb = *nodes[_id].bb.Merge(nodes[nodeId].bb);

            nodes[parB].left = _id;
            nodes[parB].right = nodeId;
            nodes[parB].parent = parA;

            nodes[_id].parent = parB;
            nodes[nodeId].parent = parB;

            if (parA == NO_USE) {
                rootId = parB;
            } else if (nodes[parA].left == _id) {
                nodes[parA].left = parB;
            } else {
                nodes[parA].right = parB;
            }

            Update(parA);
        }
    }

    std::vector<std::shared_ptr<Obj>> Search (std::shared_ptr<Obj> obj) {
        std::vector<std::shared_ptr<Obj>> found;

        std::deque<int> stack{rootId};

        const BB &bbA = obj->GetBB();

        while (!stack.empty()) {
            int s = stack.front();
            stack.pop_front();

            if (s == NO_USE) {
                continue;
            }

            const BB &bbB = nodes[s].bb;

            if (bbA.Intersects(bbB)) {
                if (nodes[s].left == NO_USE) {
                    found.emplace_back(nodes[s].obj);
                } else {
                    stack.push_back(nodes[s].left);
                    stack.push_back(nodes[s].right);
                }
            }
        }

        return found;
    }

private:
    std::vector<Node> nodes;

    int rootId;

    int CreateNode (std::shared_ptr<Obj> obj = nullptr) {
        nodes.emplace_back(obj);
        return nodes.size()-1;
    }

    void Update (int id) {
        while (id != NO_USE) {
            Node &node = nodes[id];

            node.bb = *nodes[node.left].bb.Merge(nodes[node.right].bb);

            id = node.parent;
        }
    }
};

#endif
