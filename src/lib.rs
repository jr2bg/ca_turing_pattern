use bevy::prelude::*;
/// Work from https://biologicalmodeling.org/prologue/diffusion_automaton
use rand::thread_rng;
use rand::Rng;
use rand::seq::SliceRandom;

/// CellState
/// Pair of values representing the A and B concentrations 
/// A, B in interval [0,1]
#[derive(Debug, Clone, Copy)]
pub struct CellState {
    pub a: f32,
    pub b: f32,
}

impl Default for CellState{
    fn default() -> Self {
        CellState { a: 0., b: 0. }
    }
}

impl CellState {
    fn new(a: f32, b: f32) -> Self {
        CellState { a, b }
    }

    fn change_state(&mut self, a: f32, b: f32) {
        self.a = a;
        self.b = b;
    }

    fn color(&self) -> f32 {
        if self.a + self.b <= 0. {
            return 0.
        }
        self.b / (self.a + self.b)
    }
}

/// States
/// Stores the pair of current and previous state
#[derive(Debug, Clone, Copy, Component)]
pub struct States {
    pub prev: CellState,
    pub curr: CellState,
}

impl States {
    fn initialize(a: f32, b: f32) -> Self {
        States {
            prev: CellState::default(),
            curr: CellState::new(a, b)
        }
    }

    fn shift(&mut self) {
        self.prev = self.curr;
    }

    fn get_sprite_color(&self) -> f32 {
        self.curr.color()
    }
}

/// Position
/// Pair of values indicating the row,col position of a cell
#[derive(Component)]
pub struct Position {
    pub x: usize,
    pub y: usize,
}

impl Position {
    fn new (x: usize, y: usize) -> Self {
        Position {x, y}
    }
}

/// Parameters for the simulation
/// Parameters required for the simulation of a CA for Turing patterns
/// Components:
/// `d_a` -> diffusion rate for element A in interval [0,1]
/// `d_b` -> diffusion rate for element B in interval [0,1]
/// `f` -> constant feed rate for element A in interval [0,1]
/// `k` -> constant death reaction rate for element B in interval [0,1]
/// `r` -> constant reproduction reaction rate in interval [0,1]
/// `n_rows` -> number of rows in the simulation
/// `n_cols`-> number of columns in the simulation
#[derive(Resource)]
pub struct Parameters {
    pub d_a: f32,
    pub d_b: f32,
    pub f: f32,
    pub k: f32,
    pub r: f32,
    pub n_x: usize,
    pub n_y: usize,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            d_a: 0.6,
            d_b: 0.3,
            f: 0.2,
            k: 0.1,
            r: 0.5,
            n_x: 600,
            n_y: 600,
        }
    }
}

/// Universe to be considered
/// Area where the simulation will be run
pub type Universe = Vec<Vec<CellState>>;

/// Color map
/// Area with colors for each cell
pub type ColoredMap = Vec<Vec<f32>>;

/// Initialize universe
/// Create a universe with given dimensions and n cells with
/// A and B components
pub fn initialize_universe(
    mut commands: Commands,
    parameters: Res<Parameters>,
){
    let prob_components:f64 = 0.005;

    let sprite_sz:f32 = 1.;
    let (n_x, n_y) = (parameters.n_x, parameters.n_y);
    let mut rng = rand::thread_rng();

    commands.spawn(Camera2dBundle::default());
    commands.spawn(SpatialBundle::from_transform(
        Transform::from_xyz(
            -(n_x as f32) * sprite_sz / 2., 
            -(n_y as f32) * sprite_sz / 2., 
            0.)
        ))
        .with_children(|builder| {
            for x in 0..n_x {
                for y in 0..n_y {
                    let (a, b) = if rng.gen_bool(prob_components)  {
                        (1.0_f32, 1.0_f32)
                    } else {
                        (0.0_f32, 0.0_f32)
                    };
                    builder.spawn((
                        // All elements in tuple must derive from Component
                    SpriteBundle {
                        sprite: Sprite {
                            custom_size: Some(Vec2::splat(sprite_sz)),
                            color: Color::GRAY,
                            ..Default::default()
                        },
                        transform: Transform::from_xyz(
                            sprite_sz * x as f32 ,
                            sprite_sz * y as f32,
                            0.0),
                        ..Default::default()
                    },
                    Position::new(x, y),
                    States::initialize(a, b),
                    ));
                }
            }
        });
}

/// Diffusion between two adjacent cells
/// Substract from the diffused cell the quantity of components A and B proportional to
/// its angular relation, i.e. if it is diagonal 0.05 and 0.2 in cc
/// Similar, add the corresponding quantities of A and B from the CellState at 
/// neighbour_position in  universe
fn get_adjacent_cells_diffusion(
    d_a: f32,
    d_b: f32,
    angular_rate: f32,
    diffused_cell: &mut CellState, 
    neighbour_position: Position,
    universe: &Universe
    ){

    diffused_cell.a -= angular_rate * d_a * diffused_cell.a;
    diffused_cell.b -= angular_rate * d_b * diffused_cell.b;

    diffused_cell.a += angular_rate * d_a * universe[neighbour_position.y][neighbour_position.x ].a;
    diffused_cell.b += angular_rate * d_b * universe[neighbour_position.y][neighbour_position.x ].b;
}

/// Diffusion for a cell 
/// Add the adjacent and diagonal values of substance receved due to diffusion
/// from substance A and B from its neighbours, and also substract the substance
/// given to its neighbours using `d_a` and `d_b`.
/// In this case, 0.2 and 0.05 is considered for adjacent and diagonal 
/// cells, respectively
fn get_diffusion_in_cell(
    d_a: f32,
    d_b: f32,
    cell: &CellState, 
    position: &Position,
    dimensions: &Position,
    universe: &Universe) -> CellState {

    let mut diffused_cell = *cell;

    if position.y as i32 - 1 >= 0 && position.x  as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {y: position.y - 1, x: position.x  - 1},
            universe
            );
    }

    if position.y as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {y: position.y - 1, x: position.x  },
            universe
            );
    } 

    if position.y as i32 - 1 >= 0 && position.x  + 1 < dimensions.x {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {y: position.y - 1, x: position.x  + 1},
            universe
            );
    }

    if position.x  + 1 < dimensions.x {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {y: position.y , x: position.x  + 1},
            universe
            );
    }

    if position.y + 1 < dimensions.y && position.x  + 1 < dimensions.x {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {y: position.y + 1, x: position.x  + 1},
            universe
            );
    }

    if position.y + 1 < dimensions.y {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {y: position.y + 1, x: position.x },
            universe,
            );
    }

    if position.y + 1 < dimensions.y && position.x  as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.05,
            &mut diffused_cell,
            Position {y: position.y + 1, x: position.x  - 1},
            universe
            );
    }

    if position.x  as i32 - 1 >= 0 {
        get_adjacent_cells_diffusion(
            d_a,
            d_b,
            0.2,
            &mut diffused_cell,
            Position {y: position.y , x: position.x  - 1},
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
    cell: &CellState, 
    position: &Position,
    dimensions: &Position,
    universe: &Universe,
    colored_map: &mut ColoredMap) -> CellState {

    let mut evolved_cell: CellState;

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
    
    colored_map[position.y][position.x] = color_cell(&evolved_cell);

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
    let mut evolved_universe: Universe = vec![vec![ CellState {a: 0.0, b: 0.0} ; dimensions.x]; dimensions.y];
    
    for r in 0..dimensions.y {
        for c in 0..dimensions.x{
            evolved_universe[r][c] = transition(
                parameters,
                &universe[r][c],
                &Position {y: r, x: c},
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
pub fn color_cell(cell: &CellState) -> f32 {
    if cell.a + cell.b <= 0. {
        return 0.;
    }
    cell.b / (cell.a + cell.b)
}
