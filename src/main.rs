use bevy::prelude::*;

use ca_turing_pattern::*;

fn main() {   
    // total_simulation(700, &parameters, &dimensions, universe, &mut colored_map);

    App::new()
        .insert_resource(Parameters::default())
        .insert_resource(ClearColor(Color::LIME_GREEN))
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            window: WindowDescriptor {
                title: "Turing Pattern uwu".to_string(),
                width: 900.,
                height: 700.,
                ..Default::default()
            },
            ..Default::default()
        }))
        .add_startup_system(initialize_universe)
        .run();
}
