struct Ray {
    // All elements for P and S wave
    // Distance from Earth center (km)
    double radius;
    // Direction of ray from vector pointing to Earht center (deg)
    double angle;
    // angle in radians
    double angle_rad;
    // Velocity of layer (km/s)
    double vel;
    // Spherically symmetric ray parameter (Eqn. 4.40)
    double p_sph;
    // Accumulated time of ray
    double time;
    // Distance traveled by ray (degrees)
    double dist;
};
