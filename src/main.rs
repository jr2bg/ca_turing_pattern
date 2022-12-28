use ca_turing_pattern::*;

fn main() {

    let parameters: Parameters = Parameters {
        d_a: 0.6,
        d_b: 0.3,
        f: 0.2,
        k: 0.1,
        r: 0.5,
    };
    let dimensions: Position = Position{row: 600, col: 600};
    
    let (universe , mut colored_map) = initialize_universe(&dimensions);
    
    total_simulation(700, &parameters, &dimensions, universe, &mut colored_map);

}
