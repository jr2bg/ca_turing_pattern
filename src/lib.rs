/// Work from https://biologicalmodeling.org/prologue/diffusion_automaton
use rand::thread_rng;
use rand::seq::SliceRandom;

/// Cell
/// Pair of values representing the A and B concentrations 
/// A, B in interval [0,1]
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub a: f32,
    pub b: f32,
}

/// Position
/// Pair of values indicating the row,col position of a cell
pub struct Position {
    pub row: usize,
    pub col: usize,
}

/// Universe to be considered
/// Area where the simulation will be run
pub type Universe = Vec<Vec<Cell>>;

/// Color map
/// Area with colors for each cell
pub type ColoredMap = Vec<Vec<f32>>;

/// Initialize universe
/// Create a universe with given dimensions and some values for the 
/// A and B components
pub fn initialize_universe(dimensions: &Position) -> (Universe, ColoredMap) {
    let n = 3;

    let mut universe: Universe = vec![vec![Cell {a: 0.0, b: 0.0}; dimensions.col]; dimensions.row];
    let mut colored_map: ColoredMap = vec![vec![0.0; dimensions.col]; dimensions.row];

    let mut positions: Vec<Position> = Vec::with_capacity(dimensions.row * dimensions.col);
    for r in 0..dimensions.row {
        for c in 0..dimensions.col {
            positions.push(Position{row: r, col: c});
        }
    }
    
    positions.shuffle(&mut thread_rng());
    let mut cell: &mut Cell;
    for i in 0..n {
        cell = &mut universe[positions[i].row][positions[i].col];
        *cell = Cell {a: 1.0, b: 1.0};
        colored_map[positions[i].row][positions[i].col] = color_cell(cell);
    }

    return (universe, colored_map);

}

/// Parameters for the simulation
/// Parameters required for the simulation of a CA for Turing patterns
/// Components:
/// `d_a` -> diffusion rate for element A in interval [0,1]
/// `d_b` -> diffusion rate for element B in interval [0,1]
/// `f` -> constant feed rate for element A in interval [0,1]
/// `k` -> constant death reaction rate for element B in interval [0,1]
/// `r` -> constant reproduction reaction rate
pub struct Parameters {
    pub d_a: f32,
    pub d_b: f32,
    pub f: f32,
    pub k: f32,
    pub r: f32,
}

/// Diffusion between two adjacent cells
/// Substract from the diffused cell the quantity of components A and B proportional to
/// its angular relation, i.e. if it is diagonal 0.05 and 0.2 in cc
/// Similar, add the corresponding quantities of A and B from the Cell at 
/// neighbour_position in  universe
fn get_adjacent_cells_diffusion(
    d_a: f32,
    d_b: f32,
    angular_rate: f32,
    diffused_cell: &mut Cell, 
    neighbour_position: Position,
    universe: &Universe
    ){

    diffused_cell.a -= angular_rate * d_a * diffused_cell.a;
    diffused_cell.b -= angular_rate * d_b * diffused_cell.b;

    diffused_cell.a += angular_rate * d_a * universe[neighbour_position.row][neighbour_position.col].a;
    diffused_cell.b += angular_rate * d_b * universe[neighbour_position.row][neighbour_position.col].b;
}

/// Diffusion function for each cell 
/// Add the adjacent and diagonal values of substance receved due to diffusion
/// from substance A and B from its neighbours, and also substract the substance
/// given to its neighbours using `d_a` and `d_b`.
/// In this case, 0.2 and 0.05 is considered for adjacent and diagonal 
/// cells, respectively
fn get_diffusion_in_cell(
    d_a: f32,
    d_b: f32,
    cell: &Cell, 
    position: &Position,
    dimensions: &Position,
    universe: &Universe) -> Cell {

    let mut diffused_cell = *cell;

    if position.row as i32 - 1 >= 0 && position.col as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {row: position.row - 1, col: position.col - 1},
            universe
            );
    }

    if position.row as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {row: position.row - 1, col: position.col },
            universe
            );
    } 

    if position.row as i32 - 1 >= 0 && position.col + 1 < dimensions.col {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {row: position.row - 1, col: position.col + 1},
            universe
            );
    }

    if position.col + 1 < dimensions.col {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {row: position.row, col: position.col + 1},
            universe
            );
    }

    if position.row + 1 < dimensions.row && position.col + 1 < dimensions.col {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {row: position.row + 1, col: position.col + 1},
            universe
            );
    }

    if position.row + 1 < dimensions.row {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {row: position.row + 1, col: position.col},
            universe,
            );
    }

    if position.row + 1 < dimensions.row && position.col as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {row: position.row + 1, col: position.col - 1},
            universe
            );
    }

    if position.col as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {row: position.row, col: position.col - 1},
            universe
            );
    }

    diffused_cell
}

/// Transition function
/// Considers the difussion for each cell,
/// the feed of A,
/// the death of B, and
/// the reproduction A + 2B -> 3B
fn transition(
    parameters: &Parameters,
    cell: &Cell, 
    position: &Position,
    dimensions: &Position,
    universe: &Universe,
    colored_map: &mut ColoredMap) -> Cell {

    let mut evolved_cell: Cell;

    evolved_cell = get_diffusion_in_cell(
                        parameters.d_a,
                        parameters.d_b,
                        cell,
                        position,
                        dimensions,
                        universe);

    evolved_cell.a += parameters.f * (1.0 - cell.a);

    evolved_cell.b -= parameters.k * cell.b;
    
    let reproduction_reaction: f32 = parameters.r * cell.a * cell.b.powf(2.0);
    evolved_cell.a -= reproduction_reaction;
    evolved_cell.b += reproduction_reaction;
    
    colored_map[position.row][position.col] = color_cell(&evolved_cell);

    evolved_cell
}

/// Iterate over all cells in the universe
/// From the initial state, generate another universe and return it with the
/// corresponding values of one evolution
fn evolution_universe(
    parameters: &Parameters, 
    dimensions: &Position, 
    universe: Universe,
    colored_map: &mut ColoredMap) -> Universe {
    let mut evolved_universe: Universe = vec![vec![ Cell {a: 0.0, b: 0.0} ; dimensions.col]; dimensions.row];
    
    for r in 0..dimensions.row {
        for c in 0..dimensions.col{
            evolved_universe[r][c] = transition(
                parameters,
                &universe[r][c],
                &Position {row: r, col: c},
                dimensions,
                &universe,
                colored_map
                );
        }
    }
    
    // println!("{:#?}", colored_map);

    evolved_universe
}

/// Grouped method for n-steps evolution
/// From an initial configuration of the universe, generate all the evolutions according to a given
/// n, the number of evolutions
pub fn total_simulation(
    n: i32, 
    parameters: &Parameters, 
    dimensions: &Position, 
    mut universe: Universe,
    colored_map: &mut ColoredMap){
    for _ in 0..n {
        universe = evolution_universe(
            parameters,
            dimensions,
            universe,
            colored_map
            );
    }
}

/// Color visualisation for cell
/// Give a color for each cell according to the concentrations A and B
pub fn color_cell(cell: &Cell) -> f32 {
    if cell.a + cell.b <= 0. {
        return 0.;
    }
    cell.b / (cell.a + cell.b)
}
