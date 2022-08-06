namespace Sphere
{
    internal struct Vec2D
    {
        public double m1;
        public double m2;

        public Vec2D()
        {
            m1 = 0.0f;
            m2 = 0.0f;
        }

        public double GetPolarRadius()
        {
            return Math.Sqrt(m1 * m1 + m2 * m2);
        }

        public double GetPolarAngle()
        {
            return Math.Atan2(m2, m1);
        }

        // m1 = radius : m2 = angle
        public Vec2D ToPolar()
        {
            return new Vec2D { m1 = GetPolarRadius(), m2 = GetPolarAngle() };
        }

        public Vec2D ToCart()
        {
            return new Vec2D { m1 = m1 * Math.Cos(m2), m2 = m1 * Math.Sin(m2) };
        }

        public static Vec2D operator +(Vec2D v1, Vec2D v2) => new() { m1 = v1.m1 + v2.m1, m2 = v1.m2 + v2.m2 };
        public static Vec2D operator -(Vec2D v1, Vec2D v2) => new() { m1 = v1.m1 - v2.m1, m2 = v1.m2 - v2.m2 };

        public static Vec2D operator *(Vec2D v1, double v2) => new() { m1 = v1.m1 * v2, m2 = v1.m2 * v2 };
    }

    // is actually a 2d circle
    internal class Sphere
    {
        public Vec2D pos;
        public Vec2D vel; // is in polar coordinates
        public double radius = 1.0;
        public double mass = 1.0;

        public Sphere(double x, double y)
        {
            pos = new Vec2D{ m1 = x, m2 = y};
            vel = new Vec2D{ m1 = 0.0, m2 = 0.0};
        }

        // only handles the position.
        // i recommend that when creating an array of spheres, call UpdateVel() on them individually in a triangular
        // pattern, then go across and call UpdatePos()
        public void UpdatePos(double dt)
        {
            pos += vel.ToCart() * dt;
        }

        // i recommend that when creating an array of spheres, call UpdateVel() on them individually in a triangular
        // pattern, then go across and call AdvanceSphere()
        public void UpdateVelocity(Sphere other_sphere, double dt)
        {
            Reflect(other_sphere, dt);
        }

        public double GetAnglePerpendicular(Sphere other_sphere)
        {
            // the pi / 2 part is for a 90 degree rotation
            return Math.Atan2((other_sphere.pos.m2 - pos.m2), (other_sphere.pos.m1 - pos.m1)) - (Math.PI / 2.0);
        }

        public double GetVelocityAngle()
        {
            return vel.m2;
        }

        public void ReflectOffSurface(double surface_angle)
        {
            vel.m2 = 2 * surface_angle - vel.m2;
        }

        // handles collision detection of other sphere
        public void Reflect(Sphere other_sphere, double dt)
        {
            double x_diff = (pos.m1 - other_sphere.pos.m1);
            double y_diff = (pos.m2 - other_sphere.pos.m2);
            double radius_tot = (radius + other_sphere.radius);

            if (x_diff * x_diff + y_diff * y_diff < radius_tot * radius_tot)
            {
                // ensure the spheres aren't inside of each other after colliding.
                pos -= vel.ToCart() * dt;

                ReflectOffSurface(GetAnglePerpendicular(other_sphere));
                
                // make sure to affect other sphere too.
                other_sphere.ReflectSecondary(this, GetKineticEnergy());
            }
        }

        public void ReflectSecondary(Sphere other_sphere, double other_kinetic_energy)
        {
            ReflectOffSurface(GetAnglePerpendicular(other_sphere));
            double cur_energy = GetKineticEnergy();
            
            // an elastic collision
            other_sphere.SetSpeedWithKineticEnergy(cur_energy);
            SetSpeedWithKineticEnergy(other_kinetic_energy);
        }

        // checks if the cirlce is within the circular boundary
        // true = outside boundary
        // false = inside boundary
        public bool OutsideCircularBoundary(double boundary_radius) {
            if (boundary_radius > Math.Abs(pos.GetPolarRadius()) + radius) return false;
            return true;
        }

        // handles reflection off the boundary
        // also puts the sphere inside the boundary to prevent spheres from getting stuck
        public void BoundaryReflect(double boundary_radius, double dt)
        {
            if(OutsideCircularBoundary(boundary_radius))
            {
                // gets the angle of the surface the circle is reflecting off of (the line tangent to the border circle
                // at the point at which the current circle is intersecting the border)
                double origin_angle_perpendicular = vel.m2 - Math.PI / 2.0;

                pos -= vel.ToCart() * dt; // ensure sphere within boundary

                ReflectOffSurface(origin_angle_perpendicular);

            }
        }

        public double GetKineticEnergy()
        {
            Vec2D cart = vel.ToCart();
            return 0.5 * mass * (cart.m1 * cart.m1 + cart.m2 * cart.m2);
        }

        public void SetSpeedWithKineticEnergy(double new_kinetic_energy)
        {
            vel.m1 = Math.Sqrt(2.0 * new_kinetic_energy / mass);
        }
    }

    internal class Program
    {
        static double TotalKineticEnergy(Sphere[] s)
        {
            double kinetic_energy = 0;
            foreach (Sphere sphere in s)
            {
                kinetic_energy += sphere.GetKineticEnergy();
            }

            return kinetic_energy;
        }

        static bool IsWithinCircle(Sphere[] s, double circle_radius)
        {
            foreach (Sphere sphere in s)
            {
                if (sphere.pos.GetPolarRadius() > circle_radius) return false;
            }
            return true;
        }
        
        static void Main()
        {
            File.Delete("./pos.txt");
            File.Delete("./entropy.txt");
            var rand = new Random();

            // radius of collision boundary (circle)
            const double boundary_radius = 20.0;

            // timestep of sim
            const double dt = 0.01;

            // max time spent on each simulation (recommended = 200.0, as the average simulation seems to approach a value of 100.0)
            const double max_time = 200.0;

            // the time the program should start checking the entropy (recommended = 10.0)
            const double min_check_time = 10.0;

            // debug option (false = release mode)
            const bool print_positions = false;

            // debug option (true = release mode)
            const bool do_full_sim = true;

            // the starting radius of the circle that will be used to detect the entropy (roughly)
            // if every circle is within the scan radius, then the entropy is probably lower
            double entropy_scan_radius = 20.0;

            // the rate at which to decrease the entropy radius. smaller = more continous looking graphs but also longer compute times.
            const double entropy_scan_radius_step = 0.25;

            // how many times to rerun each colliding sim to get more accurate average times for the entropy to decrease
            const int max_entropy_counter = 500;

            // list of previous times for the entropy to decrease to the radius for averaging
            List<double> previous_entropy_times = new();

            // main loop controlling the repetition of the scanning
            while (entropy_scan_radius > 0.0)
            {
                Sphere[] spheres = new Sphere[8];
                // total simulation time (for this individual sim, not total program time)
                double total_time = 0.0;

                // used for debugging
#pragma warning disable CS0219
                double last_print_time = 0.0;
#pragma warning restore CS0219

                // to spawn spheres, we will use the equation r = 2.2 * theta to find where to put the center of the circle
                // with theta increasing by 1 for the next position of the circle.
                Vec2D polar_circle_point = new() { m1 = 0, m2 = 0 };

                // sphere initialization
                for (int i = 0; i < spheres.Length; i++)
                {
                    Vec2D sphere_center = polar_circle_point.ToCart();
                    spheres[i] = new Sphere(sphere_center.m1, sphere_center.m2)
                    {
                        // give them random velocities
                        // this doesn't choose a point over a circle with uniform probability, but it doesn't need to. It's good enough
                        vel = new Vec2D { m1 = rand.NextDouble() * 12 - 6, m2 = rand.NextDouble() * 12 - 6 }.ToPolar()
                    };

                    polar_circle_point.m2 += 1.0; // theta

                    // if the circle would be placed outside of the boundary
                    if (2.2 * polar_circle_point.m2 > boundary_radius - 1.0)
                    {
                        throw new Exception("ERROR: boundary radius too small, or too many circles");
                    }
                    polar_circle_point.m1 = 2.2 * polar_circle_point.m2;
                }

                // the loop for the individual simulations
                while (total_time < max_time)
                {
                    // debug code to print position of individual circle
                    if (print_positions && (last_print_time == 0.0 || total_time - last_print_time > 0.5))
                    {
                        File.AppendAllText("./pos.txt", total_time.ToString("F3") + '\t' + (spheres[0].pos.m1).ToString("F3") + "\t" + TotalKineticEnergy(spheres).ToString("F3") + "\n");
                        last_print_time = total_time;
                    }

                    // updating velocity, collision detection
                    for (int j = 0; j < spheres.Length - 1; j++)
                    {
                        // updates every circle boundary detection except for the last circle in the list.
                        spheres[j].BoundaryReflect(boundary_radius, dt);
                        for (int k = j + 1; k < spheres.Length; k++)
                        {
                            // we use mod to make sure the first circle is favored too much by reversing directions every other iteration
                            // this probably doesn't matter very much
                            if (j % 2 == 0)
                            {
                                spheres[j].UpdateVelocity(spheres[^k], dt);
                            }
                            else
                            {
                                spheres[j].UpdateVelocity(spheres[k], dt);
                            }
                        }
                    }

                    // the loop above misses the last sphere
                    spheres[^1].BoundaryReflect(boundary_radius, dt);

                    // updating position
                    for (int j = 0; j < spheres.Length; j++)
                    {
                        spheres[j].UpdatePos(dt);
                    }

                    // check the entropy with the circle method
                    // the total_time check is to allow the system settle a bit.
                    if (do_full_sim && total_time > min_check_time)
                    {
                        if (IsWithinCircle(spheres, entropy_scan_radius))
                        {
                            if (previous_entropy_times.Count > max_entropy_counter)
                            {
                                double average = previous_entropy_times.Sum();
                                average /= previous_entropy_times.Count;
                                previous_entropy_times.Clear();

                                entropy_scan_radius -= entropy_scan_radius_step;
                                Console.WriteLine("Current radius difference: " + (boundary_radius - entropy_scan_radius).ToString("F3"));

                                File.AppendAllText("./entropy.txt", (boundary_radius - entropy_scan_radius).ToString("F4") + '\t' + average.ToString("F4") + "\n");
                            }
                            else
                            {
                                previous_entropy_times.Add(total_time);
                            }
                            // make sure to move on to next iteration, otherwise the program will keep simulating the current simulation.
                            break;
                        }
                    }

                    total_time += dt;
                }

                // enabled in debug situations
#pragma warning disable CS0162 
                if (!do_full_sim) goto exit_loop;
#pragma warning restore CS0162 // Unreachable code detected
            }
        exit_loop:;
        }
        
    }
}