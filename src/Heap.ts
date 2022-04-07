// https://github.com/qiao/heap.js

/*
Default comparison to be used
 */
function defaultCmp<T>(x: T, y: T): number {
    if (x < y) { return -1; }
    if (x > y) { return 1; }
    return 0;
}

export default class Heap<T> {

    /*
    Insert item x in list a, and keep it sorted assuming a is sorted.
    
    If x is already in a, insert it to the right of the rightmost x.
    
    Optional args lo (default 0) and hi (default a.length) bound the slice
    of a to be searched.
     */
    private insort<T>(array: T[], x: T, low: number, high: number | null, cmp: (x: T, y: T) => number = defaultCmp): void {

        if (low == null) { low = 0; }
        if (low < 0) { throw new Error('lo must be non-negative'); }

        if (high == null) { high = array.length; }

        let mid: number;
        while (low < high) {
            mid = Math.floor((low + high) / 2);
            if (cmp(x, array[mid]) < 0) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }

        array.splice(low, 0, x);
    }

    /*
    Push item onto heap, maintaining the heap inletiant.
     */
    private _heappush<T>(array: T[], item: T, cmp: (x: T, y: T) => number = defaultCmp): void {
        array.push(item);
        this._siftdown(array, 0, array.length - 1, cmp);
    }

    /*
    Pop the smallest item off the heap, maintaining the heap inletiant.
     */
    private _heappop<T>(array: T[], cmp: (x: T, y: T) => number = defaultCmp): T | null {
        if (!array.length) return null;

        let lastelt: T = array.pop() as T;
        let returnitem: T;

        if (array.length) {
            returnitem = array[0];
            array[0] = lastelt;
            this._siftup(array, 0, cmp);
        } else {
            returnitem = lastelt;
        }
        return returnitem;
    }

    /*
    Pop and return the current smallest value, and add the new item.
    
    This is more efficient than heappop() followed by heappush(), and can be
    more appropriate when using a fixed size heap. Note that the value
    returned may be larger than item! That constrains reasonable use of
    this routine unless written as part of a conditional replacement:
        if item > array[0]
          item = heapreplace(array, item)
     */
    private heapreplace<T>(array: T[], item: T, cmp: (x: T, y: T) => number = defaultCmp): T {
        let returnitem: T = array[0];
        array[0] = item;
        this._siftup(array, 0, cmp);
        return returnitem;
    }

    /*
    Fast version of a heappush followed by a heappop.
     */
    private heappushpop<T>(array: T[], item: T, cmp: (x: T, y: T) => number = defaultCmp): T {
        if (array.length && cmp(array[0], item) < 0) {
            array[0] = item;
            this._siftup(array, 0, cmp);
        }
        return item;
    }

    /*
    Transform list into a heap, in-place, in O(array.length) time.
     */
    private _heapify<T>(array: T[], cmp: (x: T, y: T) => number = defaultCmp): void {
        for (let i = Math.floor(array.length / 2) - 1; i >= 0; i--) {
            this._siftup(array, i, cmp)
        }
    }

    /*
    Update the position of the given item in the heap.
    This should be called every time the item is being modified.
     */
    private _updateItem<T>(array: T[], item: T, cmp: (x: T, y: T) => number = defaultCmp): void {
        let pos = array.indexOf(item);
        if (pos === -1) return;
        this._siftdown(array, 0, pos, cmp);
        this._siftup(array, pos, cmp);
    }

    /*
    Find the n largest elements in a dataset.
     */
    private nlargest<T>(array: T[], n: number, cmp: (x: T, y: T) => number = defaultCmp): T[] {
        let result: T[] = array.slice(0, n);
        if (!result.length) {
            return result;
        }
        this._heapify(result, cmp);

        for (let i = n; i < array.length; i++) {
            this.heappushpop(result, array[i], cmp);
        }
        return result.sort(cmp).reverse();
    }

    /*
    Find the n smallest elements in a dataset.
     */

    private nsmallest<T>(array: T[], n: number, cmp: (x: T, y: T) => number = defaultCmp): T[] {
        let result: T[];
        if (n * 10 <= array.length) {
            result = array.slice(0, n).sort(cmp);
            if (!result.length) return [];

            let los = result[result.length - 1];
            for (let i = n; i < array.length; i++) {
                if (cmp(array[i], los) < 0) {
                    this.insort(result, array[i], 0, null, cmp);
                    result.pop();
                    los = result[result.length - 1];
                }
            }
            return result;
        }

        this._heapify(array, cmp);

        result = [];
        for (let i = 0; i < Math.min(n, array.length); i++) {
            result.push(this._heappop(array, cmp) as T);
        }
        return result;
    }

    private _siftdown<T>(array: T[], startpos: number, pos: number, cmp: (x: T, y: T) => number = defaultCmp): void {
        let newitem = array[pos];
        let parent, parentpos;
        while (pos > startpos) {
            parentpos = (pos - 1) >> 1;
            parent = array[parentpos];
            if (cmp(newitem, parent) < 0) {
                array[pos] = parent;
                pos = parentpos;
                continue;
            }
            break;
        }
        array[pos] = newitem;
    }

    private _siftup<T>(array: T[], pos: number, cmp: (x: T, y: T) => number = defaultCmp): void {
        let childpos, endpos, newitem, rightpos, startpos;
        if (cmp == null) {
            cmp = defaultCmp;
        }
        endpos = array.length;
        startpos = pos;
        newitem = array[pos];
        childpos = 2 * pos + 1;
        while (childpos < endpos) {
            rightpos = childpos + 1;
            if (rightpos < endpos && !(cmp(array[childpos], array[rightpos]) < 0)) {
                childpos = rightpos;
            }
            array[pos] = array[childpos];
            pos = childpos;
            childpos = 2 * pos + 1;
        }
        array[pos] = newitem;
        return this._siftdown(array, startpos, pos, cmp);
    }


    private _nodes: T[] = [];

    private _cmp: (x: T, y: T) => number;

    constructor(cmp: (x: T, y: T) => number = defaultCmp) {
        this._cmp = cmp;
        this._nodes = [];
    }

    public push(x: T) {
        this._heappush(this._nodes, x, this._cmp);
    }

    public pop(): T | null {
        return this._heappop(this._nodes, this._cmp);
    }

    public peek(): T {
        return this._nodes[0];
    }

    public contains(x: T): boolean {
        return this._nodes.indexOf(x) !== -1;
    }

    public replace(x: T): T {
        return this.heapreplace(this._nodes, x, this._cmp);
    }

    public pushpop(x: T): T {
        return this.heappushpop(this._nodes, x, this._cmp);
    }

    public heapify(): void {
        this._heapify(this._nodes, this._cmp);
    }

    public updateItem(x: T): void {
        this._updateItem(this._nodes, x, this._cmp);
    }

    public clear(): void {
        this._nodes = [];
    }

    public get empty(): boolean {
        return this._nodes.length === 0;
    }

    public get size(): number {
        return this._nodes.length;
    }

    public clone(): Heap<T> {
        let heap: Heap<T> = new Heap();
        heap._nodes = this._nodes.slice(0);
        return heap;
    }

    public toArray(): T[] {
        return this._nodes.slice(0);
    }
}

