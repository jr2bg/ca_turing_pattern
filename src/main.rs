use bevy::prelude::*;

use ca_turing_pattern::*;

fn main() {

    let parameters: Parameters = Parameters {
        d_a: 0.6,
        d_b: 0.3,
        f: 0.2,
        k: 0.1,
        r: 0.5,
        n_rows: 600,
        n_cols: 600,
    };
    let dimensions: Position = Position{row: 600, col: 600};
    
    let (universe , mut colored_map) = initialize_universe(&dimensions);
    
    total_simulation(700, &parameters, &dimensions, universe, &mut colored_map);

    App::new()
        .insert_resource(parameters)
        .insert_resource(ClearColor(Color::LIME_GREEN))
        .add_plugins(DefaultPlugins)
        .add_startup_system(initialize_universe);
}
