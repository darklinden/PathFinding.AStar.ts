import PathTestScenarios from "./PathTestScenarios";
const scenarios = PathTestScenarios;

import pathfinding from './PathFinding';
type AStar = pathfinding.AStar;

function test(finder: AStar, startX: number, startY: number, endX: number, endY: number, grid: pathfinding.Grid) {
    let path = finder.findPath(startX, startY, endX, endY, grid);
    console.log(grid.toColorString([startX, startY], [endX, endY], path));
    console.log(path);
    console.log('path[0] == [startX, startY]', path[0] == [startX, startY]);
    console.log('path[path.length - 1] == [endX, endY]', path[path.length - 1] == [endX, endY]);
};

function pathTest() {
    const finder = new pathfinding.AStar();
    // Load all the scenarios and test against the finder.
    for (let i = 0; i < scenarios.length; ++i) {
        let scen = scenarios[i];
        let matrix = scen.matrix;
        let height = matrix.length;
        let width = matrix[0].length;
        let grid = new pathfinding.Grid(width, height, matrix);
        test(
            finder,
            scen.startX, scen.startY,
            scen.endX, scen.endY,
            grid
        );
    }
}

pathTest();
