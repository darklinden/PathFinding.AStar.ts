import Heap from "./Heap";

namespace pathfinding {

    export enum DiagonalMovement {
        Always = 1,
        Never = 2,
        IfAtMostOneObstacle = 3,
        OnlyWhenNoObstacles = 4
    }

    const F: number = Math.SQRT2 - 1;

    export class Heuristic {

        /**
           * Manhattan distance.
           */
        public static manhattan(dx: number, dy: number): number {
            return dx + dy;
        }

        /**
         * Euclidean distance.
         */
        public static euclidean(dx: number, dy: number): number {
            return Math.sqrt(dx * dx + dy * dy);
        }

        /**
         * Octile distance.
         */
        public static octile(dx: number, dy: number): number {
            return (dx < dy) ? F * dx + dy : F * dy + dx;
        }

        /**
         * Chebyshev distance.
         */
        public static chebyshev(dx: number, dy: number): number {
            return Math.max(dx, dy);
        }
    }

    export class Node {
        public x: number;
        public y: number;
        public walkable: boolean;
        public parent: Node | null = null;
        public g: number = 0;
        public f: number = 0;
        public h: number = 0;
        public opened: boolean = false;
        public closed: boolean = false;

        constructor(x: number, y: number, walkable: boolean = true) {
            this.x = x;
            this.y = y;
            this.walkable = walkable;
        }
    }

    /**
     * The Grid class, which serves as the encapsulation of the layout of the nodes.
     */
    export class Grid {

        public width: number;
        public height: number;
        public nodes: Node[][];

        /**
         * The Grid class, which serves as the encapsulation of the layout of the nodes.
         * @constructor
         * @param {number|Array.<Array.<(number|boolean)>>} width_or_matrix Number of columns of the grid, or matrix
         * @param {number} height Number of rows of the grid.
         * @param {Array.<Array.<(number|boolean)>>} [matrix] - A 0-1 matrix
         *     representing the walkable status of the nodes(0 or false for walkable).
         *     If the matrix is not supplied, all the nodes will be walkable.
         */
        constructor(width_or_matrix: number | Array<Array<(number | boolean)>>, height: number, matrix: Array<Array<(number | boolean)>> | null) {
            let width: number;

            if (typeof width_or_matrix !== 'object') {
                width = width_or_matrix;
            } else {
                height = width_or_matrix.length;
                width = width_or_matrix[0].length;
                matrix = width_or_matrix;
            }

            this.width = width;
            this.height = height;
            this.nodes = this._buildNodes(width, height, matrix);
        }

        /**
         * Build and return the nodes.
         * @private
         * @param {number} width
         * @param {number} height
         * @param {Array.<Array.<number|boolean>>} [matrix] - A 0-1 matrix representing the walkable status of the nodes.
         * @see Grid
         */
        _buildNodes(width: number, height: number, matrix: Array<Array<number | boolean>> | null): Node[][] {
            let i, j, nodes = new Array(height);

            for (i = 0; i < height; ++i) {
                nodes[i] = new Array(width);
                for (j = 0; j < width; ++j) {
                    nodes[i][j] = new Node(j, i);
                }
            }

            if (!matrix) {
                return nodes;
            }

            if (matrix.length !== height || matrix[0].length !== width) {
                throw new Error('Matrix size does not fit');
            }

            for (i = 0; i < height; ++i) {
                for (j = 0; j < width; ++j) {
                    if (matrix[i][j]) {
                        // 0, false, null will be walkable
                        // while others will be un-walkable
                        nodes[i][j].walkable = false;
                    }
                }
            }

            return nodes;
        }

        public getNodeAt(x: number, y: number) {
            return this.nodes[y][x];
        }

        /**
         * Determine whether the node at the given position is walkable.
         * (Also returns false if the position is outside the grid.)
         */
        public isWalkableAt(x: number, y: number): boolean {
            return this.isInside(x, y) && this.nodes[y][x].walkable;
        }

        /**
         * Determine whether the position is inside the grid.
         * XXX: `grid.isInside(x, y)` is wierd to read.
         * It should be `(x, y) is inside grid`, but I failed to find a better name for this method.
         */
        public isInside(x: number, y: number): boolean {
            return (x >= 0 && x < this.width) && (y >= 0 && y < this.height);
        }

        /**
         * Set whether the node on the given position is walkable.
         * NOTE: throws exception if the coordinate is not inside the grid.
         */
        public setWalkableAt(x: number, y: number, walkable: boolean): void {
            this.nodes[y][x].walkable = walkable;
        }

        /**
         * Get the neighbors of the given node.
         *
         *     offsets      diagonalOffsets:
         *  +---+---+---+    +---+---+---+
         *  |   | 0 |   |    | 0 |   | 1 |
         *  +---+---+---+    +---+---+---+
         *  | 3 |   | 1 |    |   |   |   |
         *  +---+---+---+    +---+---+---+
         *  |   | 2 |   |    | 3 |   | 2 |
         *  +---+---+---+    +---+---+---+
         *
         *  When allowDiagonal is true, if offsets[i] is valid, then
         *  diagonalOffsets[i] and
         *  diagonalOffsets[(i + 1) % 4] is valid.
         */
        public getNeighbors(node: Node, diagonalMovement: DiagonalMovement): Node[] {
            let x = node.x, y = node.y, neighbors = [], s0 = false, d0 = false, s1 = false, d1 = false, s2 = false, d2 = false, s3 = false, d3 = false, nodes = this.nodes;

            // ↑
            if (this.isWalkableAt(x, y - 1)) {
                neighbors.push(nodes[y - 1][x]);
                s0 = true;
            }
            // →
            if (this.isWalkableAt(x + 1, y)) {
                neighbors.push(nodes[y][x + 1]);
                s1 = true;
            }
            // ↓
            if (this.isWalkableAt(x, y + 1)) {
                neighbors.push(nodes[y + 1][x]);
                s2 = true;
            }
            // ←
            if (this.isWalkableAt(x - 1, y)) {
                neighbors.push(nodes[y][x - 1]);
                s3 = true;
            }

            if (diagonalMovement === DiagonalMovement.Never) {
                return neighbors;
            }

            if (diagonalMovement === DiagonalMovement.OnlyWhenNoObstacles) {
                d0 = s3 && s0;
                d1 = s0 && s1;
                d2 = s1 && s2;
                d3 = s2 && s3;
            } else if (diagonalMovement === DiagonalMovement.IfAtMostOneObstacle) {
                d0 = s3 || s0;
                d1 = s0 || s1;
                d2 = s1 || s2;
                d3 = s2 || s3;
            } else if (diagonalMovement === DiagonalMovement.Always) {
                d0 = true;
                d1 = true;
                d2 = true;
                d3 = true;
            } else {
                throw new Error('Incorrect value of diagonalMovement');
            }

            // ↖
            if (d0 && this.isWalkableAt(x - 1, y - 1)) {
                neighbors.push(nodes[y - 1][x - 1]);
            }
            // ↗
            if (d1 && this.isWalkableAt(x + 1, y - 1)) {
                neighbors.push(nodes[y - 1][x + 1]);
            }
            // ↘
            if (d2 && this.isWalkableAt(x + 1, y + 1)) {
                neighbors.push(nodes[y + 1][x + 1]);
            }
            // ↙
            if (d3 && this.isWalkableAt(x - 1, y + 1)) {
                neighbors.push(nodes[y + 1][x - 1]);
            }

            return neighbors;
        }

        /**
         * Get a clone of this grid.
         */
        public clone(): Grid {
            const newGrid = new Grid(this.width, this.height, null);

            for (let i = 0; i < this.height; ++i) {
                for (let j = 0; j < this.width; ++j) {
                    newGrid.nodes[i][j].walkable = this.nodes[i][j].walkable;
                }
            }

            return newGrid;
        }

        /**
         * Get a string map
         */
        public toString(start: [number, number] | null = null, end: [number, number] | null = null, road: Array<[number, number]> | null = null): string {
            let map = 'map:\n';
            for (let i = 0; i < this.height; ++i) {
                for (let j = 0; j < this.width; ++j) {
                    if (start && j == start[0] && i == start[1]) {
                        map += 'S, ';
                    }
                    else if (end && j == end[0] && i == end[1]) {
                        map += 'E, ';
                    }
                    else if (road && this.nodes[i][j].walkable && road.findIndex(e => e[0] == j && e[1] == i) != -1) {
                        map += 'R, ';
                    }
                    else {
                        map += this.nodes[i][j].walkable ? '1, ' : '0, ';
                    }
                }
                map += '\n';
            }
            return map;
        }

        /**
         * Get a string map
         */
        public toColorString(start: [number, number] | null = null, end: [number, number] | null = null, road: Array<[number, number]> | null = null): string {

            const Reset = "\x1b[0m";
            const Red = "\x1b[31m";
            const Green = "\x1b[32m";
            const Yellow = "\x1b[33m";
            const Blue = "\x1b[34m";

            let map = 'map:\n';
            for (let i = 0; i < this.height; ++i) {
                for (let j = 0; j < this.width; ++j) {
                    if (start && j == start[0] && i == start[1]) {
                        map += Red + 'S, ' + Reset;
                    }
                    else if (end && j == end[0] && i == end[1]) {
                        map += Blue + 'E, ' + Reset;
                    }
                    else if (road && this.nodes[i][j].walkable && road.findIndex(e => e[0] == j && e[1] == i) != -1) {
                        map += Green + 'R, ' + Reset;
                    }
                    else {
                        map += this.nodes[i][j].walkable ? '1, ' : '0, ';
                    }
                }
                map += '\n';
            }
            return map;
        }
    }

    /**
     * Backtrace according to the parent records and return the path.
     * (including both start and end nodes)
     */
    export function backtrace(node: Node): Array<[number, number]> {
        let path: Array<[number, number]> = [[node.x, node.y]];
        while (node.parent) {
            node = node.parent;
            path.push([node.x, node.y]);
        }
        return path.reverse();
    }

    /**
     * Backtrace from start and end node, and return the path.
     * (including both start and end nodes)
     */
    export function biBacktrace(nodeA: Node, nodeB: Node) {
        var pathA = backtrace(nodeA),
            pathB = backtrace(nodeB);
        return pathA.concat(pathB.reverse());
    }

    /**
     * Compute the length of the path.
     */
    export function pathLength(path: Array<[number, number]>): number {
        var i, sum = 0, a, b, dx, dy;
        for (i = 1; i < path.length; ++i) {
            a = path[i - 1];
            b = path[i];
            dx = a[0] - b[0];
            dy = a[1] - b[1];
            sum += Math.sqrt(dx * dx + dy * dy);
        }
        return sum;
    }

    /**
     * Given the start and end coordinates, return all the coordinates lying
     * on the line formed by these coordinates, based on Bresenham's algorithm.
     * http://en.wikipedia.org/wiki/Bresenham's_line_algorithm#Simplification
     * @param {number} x0 Start x coordinate
     * @param {number} y0 Start y coordinate
     * @param {number} x1 End x coordinate
     * @param {number} y1 End y coordinate
     * @return {Array.<Array.<number>>} The coordinates on the line
     */
    export function interpolate(x0: number, y0: number, x1: number, y1: number): Array<[number, number]> {
        const line: Array<[number, number]> = [];

        let sx, sy, dx, dy, err, e2: number;

        dx = Math.abs(x1 - x0);
        dy = Math.abs(y1 - y0);

        sx = (x0 < x1) ? 1 : -1;
        sy = (y0 < y1) ? 1 : -1;

        err = dx - dy;

        while (true) {
            line.push([x0, y0]);

            if (x0 === x1 && y0 === y1) {
                break;
            }

            e2 = 2 * err;
            if (e2 > -dy) {
                err = err - dy;
                x0 = x0 + sx;
            }
            if (e2 < dx) {
                err = err + dx;
                y0 = y0 + sy;
            }
        }

        return line;
    }

    /**
     * Given a compressed path, return a new path that has all the segments
     * in it interpolated.
     * @param {Array.<Array.<number>>} path The path
     * @return {Array.<Array.<number>>} expanded path
     */
    export function expandPath(path: Array<[number, number]>): Array<[number, number]> {
        const expanded: Array<[number, number]> = [];
        let
            len = path.length,
            coord0, coord1,
            interpolated,
            interpolatedLen,
            i, j;

        if (len < 2) {
            return expanded;
        }

        for (i = 0; i < len - 1; ++i) {
            coord0 = path[i];
            coord1 = path[i + 1];

            interpolated = interpolate(coord0[0], coord0[1], coord1[0], coord1[1]);
            interpolatedLen = interpolated.length;
            for (j = 0; j < interpolatedLen - 1; ++j) {
                expanded.push(interpolated[j]);
            }
        }
        expanded.push(path[len - 1]);

        return expanded;
    }

    /**
     * Smoothen the give path.
     * The original path will not be modified; a new path will be returned.
     * @param {PF.Grid} grid
     * @param {Array.<Array.<number>>} path The path
     */
    export function smoothenPath(grid: Grid, path: Array<[number, number]>): Array<[number, number]> {
        let len: number = path.length,
            x0: number = path[0][0],        // path start x
            y0: number = path[0][1],        // path start y
            x1: number = path[len - 1][0],  // path end x
            y1: number = path[len - 1][1],  // path end y
            sx: number, sy: number,                 // current start coordinate
            ex: number, ey: number,                 // current end coordinate
            i: number, j: number,
            coord: [number, number],
            line: Array<[number, number]>,
            testCoord: [number, number],
            blocked: boolean,
            lastValidCoord: [number, number];

        sx = x0;
        sy = y0;
        const newPath: Array<[number, number]> = [[sx, sy]];

        for (i = 2; i < len; ++i) {
            coord = path[i];
            ex = coord[0];
            ey = coord[1];
            line = interpolate(sx, sy, ex, ey);

            blocked = false;
            for (j = 1; j < line.length; ++j) {
                testCoord = line[j];

                if (!grid.isWalkableAt(testCoord[0], testCoord[1])) {
                    blocked = true;
                    break;
                }
            }
            if (blocked) {
                lastValidCoord = path[i - 1];
                newPath.push(lastValidCoord);
                sx = lastValidCoord[0];
                sy = lastValidCoord[1];
            }
        }
        newPath.push([x1, y1]);

        return newPath;
    }

    /**
     * Compress a path, remove redundant nodes without altering the shape
     * The original path is not modified
     * @param {Array.<Array.<number>>} path The path
     * @return {Array.<Array.<number>>} The compressed path
     */
    export function compressPath(path: Array<[number, number]>): Array<[number, number]> {

        // nothing to compress
        if (path.length < 3) {
            return path;
        }

        var compressed: Array<[number, number]> = [],
            sx: number = path[0][0], // start x
            sy: number = path[0][1], // start y
            px: number = path[1][0], // second point x
            py: number = path[1][1], // second point y
            dx: number = px - sx, // direction between the two points
            dy: number = py - sy, // direction between the two points
            lx: number, ly: number,
            ldx: number, ldy: number,
            sq: number, i: number;

        // normalize the direction
        sq = Math.sqrt(dx * dx + dy * dy);
        dx /= sq;
        dy /= sq;

        // start the new path
        compressed.push([sx, sy]);

        for (i = 2; i < path.length; i++) {

            // store the last point
            lx = px;
            ly = py;

            // store the last direction
            ldx = dx;
            ldy = dy;

            // next point
            px = path[i][0];
            py = path[i][1];

            // next direction
            dx = px - lx;
            dy = py - ly;

            // normalize
            sq = Math.sqrt(dx * dx + dy * dy);
            dx /= sq;
            dy /= sq;

            // if the direction has changed, store the point
            if (dx !== ldx || dy !== ldy) {
                compressed.push([lx, ly]);
            }
        }

        // store the last point
        compressed.push([px, py]);

        return compressed;
    }

    /**
     * IAStarOpt
     * @param {DiagonalMovement} opt.diagonalMovement Allowed diagonal movement.
     * @param {function} opt.heuristic Heuristic function to estimate the distance (defaults to manhattan).
     * @param {integer} opt.weight Weight to apply to the heuristic to allow for suboptimal paths,
     */
    export interface IAStarOpt {
        diagonalMovement: DiagonalMovement;
        heuristic: (dx: number, dy: number) => number;
        weight: number;
    }

    const defaultAStarOpt: IAStarOpt = {
        diagonalMovement: DiagonalMovement.Always,
        heuristic: Heuristic.manhattan,
        weight: 1
    };

    /**
     * A* path-finder.
     * based upon https://github.com/bgrins/javascript-astar
     */
    export class AStar {

        private heuristic: (dx: number, dy: number) => number;
        private weight: number = 1;
        private diagonalMovement: DiagonalMovement = DiagonalMovement.Always;

        constructor(opt: IAStarOpt = defaultAStarOpt) {
            this.heuristic = opt.heuristic || Heuristic.manhattan;
            this.weight = opt.weight || 1;
            this.diagonalMovement = opt.diagonalMovement;

            // When diagonal movement is allowed the manhattan heuristic is not admissible
            // It should be octile instead
            if (this.diagonalMovement === DiagonalMovement.Never) {
                this.heuristic = opt.heuristic || Heuristic.manhattan;
            } else {
                this.heuristic = opt.heuristic || Heuristic.octile;
            }
        }

        /**
         * Find and return the the path.
         * @return {Array.<[number, number]>} The path, including both start and
         *     end positions.
         */
        public findPath(startX: number, startY: number, endX: number, endY: number, grid: Grid): Array<[number, number]> {
            const openList = new Heap((nodeA: Node, nodeB: Node) => { return nodeA.f - nodeB.f; });
            const startNode = grid.getNodeAt(startX, startY);
            const endNode = grid.getNodeAt(endX, endY);
            let node: Node;
            let neighbors: Node[];

            //  Math.abs
            //   SQRT2 = Math.SQRT2
            //    node, neighbors, neighbor, i, l, x, y, ng;

            // set the `g` and `f` value of the start node to be 0
            startNode.g = 0;
            startNode.f = 0;

            // push the start node into the open list
            openList.push(startNode);
            startNode.opened = true;

            // while the open list is not empty
            while (!openList.empty) {
                // pop the position of node which has the minimum `f` value.
                node = openList.pop() as Node;
                node.closed = true;

                // if reached the end position, construct the path and return it
                if (node === endNode) {
                    return backtrace(endNode);
                }

                // get neigbours of the current node
                neighbors = grid.getNeighbors(node, this.diagonalMovement);
                for (let i = 0, l = neighbors.length; i < l; ++i) {
                    const neighbor = neighbors[i];

                    if (neighbor.closed) {
                        continue;
                    }

                    const x = neighbor.x;
                    const y = neighbor.y;

                    // get the distance between current node and the neighbor
                    // and calculate the next g score
                    const ng = node.g + ((x - node.x === 0 || y - node.y === 0) ? 1 : Math.SQRT2);

                    // check if the neighbor has not been inspected yet, or
                    // can be reached with smaller cost from the current node
                    if (!neighbor.opened || ng < neighbor.g) {
                        neighbor.g = ng;
                        neighbor.h = neighbor.h || this.weight * this.heuristic(Math.abs(x - endX), Math.abs(y - endY));
                        neighbor.f = neighbor.g + neighbor.h;
                        neighbor.parent = node;

                        if (!neighbor.opened) {
                            openList.push(neighbor);
                            neighbor.opened = true;
                        } else {
                            // the neighbor can be reached with smaller cost.
                            // Since its f value has been updated, we have to
                            // update its position in the open list
                            openList.updateItem(neighbor);
                        }
                    }
                } // end for each neighbor
            }     // end while not open list empty

            // fail to find the path
            return [];
        }
    }
}

export default pathfinding;