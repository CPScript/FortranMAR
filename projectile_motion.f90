program projectile_motion
    implicit none
    
    ! Physical constants and parameters
    real(8), parameter :: g = 9.81d0          ! Gravitational acceleration (m/s²)
    real(8), parameter :: rho_air = 1.225d0   ! Air density at sea level (kg/m³)
    real(8), parameter :: pi = 3.141592653589793d0
    
    ! Projectile parameters
    real(8) :: mass = 0.145d0                 ! Baseball mass (kg)
    real(8) :: radius = 0.037d0               ! Baseball radius (m)
    real(8) :: drag_coeff = 0.47d0            ! Drag coefficient for sphere
    real(8) :: cross_area                     ! Cross-sectional area
    real(8) :: drag_factor                    ! Combined drag factor
    
    ! Initial conditions
    real(8) :: v0 = 45.0d0                    ! Initial velocity (m/s)
    real(8) :: angle_deg = 35.0d0             ! Launch angle (degrees)
    real(8) :: angle_rad                      ! Launch angle (radians)
    real(8) :: x, y                           ! Position coordinates
    real(8) :: vx, vy                         ! Velocity components
    
    ! Numerical integration parameters
    real(8) :: dt = 0.001d0                   ! Time step (s)
    real(8) :: t = 0.0d0                      ! Current time
    real(8) :: max_time = 20.0d0              ! Maximum simulation time
    integer :: step_count = 0
    integer, parameter :: output_interval = 100 ! Output every N steps
    
    ! File handling
    integer, parameter :: file_unit = 10
    character(len=50) :: filename = "trajectory_data.csv"
    
    ! Calculate derived parameters
    cross_area = pi * radius**2
    drag_factor = 0.5d0 * rho_air * cross_area * drag_coeff / mass
    angle_rad = angle_deg * pi / 180.0d0
    
    ! Set initial conditions
    x = 0.0d0
    y = 0.0d0
    vx = v0 * cos(angle_rad)
    vy = v0 * sin(angle_rad)
    
    ! Open output file
    open(unit=file_unit, file=filename, status='replace', action='write')
    write(file_unit, '(A)') "time,x_pos,y_pos,velocity_x,velocity_y,total_velocity"
    
    ! Display simulation parameters
    write(*, '(A)') "Projectile Motion Simulation with Air Resistance"
    write(*, '(A, F8.3, A)') "Initial Velocity: ", v0, " m/s"
    write(*, '(A, F6.1, A)') "Launch Angle: ", angle_deg, " degrees"
    write(*, '(A, F8.6, A)') "Mass: ", mass, " kg"
    write(*, '(A, F8.6, A)') "Drag Coefficient: ", drag_coeff
    write(*, '(A, F8.3, A)') "Time Step: ", dt, " s"
    write(*, '(A)') ""
    write(*, '(A)') "Time(s)    X(m)      Y(m)      Vx(m/s)   Vy(m/s)   |V|(m/s)"
    write(*, '(A)') "-----------------------------------------------------------"
    
    ! Main simulation loop using 4th-order Runge-Kutta integration
    do while (y >= 0.0d0 .and. t <= max_time)
        
        ! Output current state
        if (mod(step_count, output_interval) == 0) then
            write(*, '(F6.3, 5F10.3)') t, x, y, vx, vy, sqrt(vx**2 + vy**2)
        end if
        
        ! Write to file
        write(file_unit, '(F12.6, A, F12.6, A, F12.6, A, F12.6, A, F12.6, A, F12.6)') &
            t, ",", x, ",", y, ",", vx, ",", vy, ",", sqrt(vx**2 + vy**2)
        
        ! Runge-Kutta 4th order integration
        call runge_kutta_step(x, y, vx, vy, dt, drag_factor)
        
        t = t + dt
        step_count = step_count + 1
    end do
    
    ! Final output
    write(*, '(A)') "--------------------------------------------------------"
    write(*, '(A, F8.3, A)') "Flight time: ", t, " seconds"
    write(*, '(A, F8.3, A)') "Range: ", x, " meters"
    write(*, '(A, F8.3, A)') "Max height reached during flight"
    write(*, '(A, A)') "Trajectory data saved to: ", filename
    
    close(file_unit)
    
contains
    
    ! Runge-Kutta 4th order integration subroutine
    subroutine runge_kutta_step(x, y, vx, vy, dt, drag_factor)
        implicit none
        real(8), intent(inout) :: x, y, vx, vy
        real(8), intent(in) :: dt, drag_factor
        
        real(8) :: k1x, k1y, k1vx, k1vy
        real(8) :: k2x, k2y, k2vx, k2vy  
        real(8) :: k3x, k3y, k3vx, k3vy
        real(8) :: k4x, k4y, k4vx, k4vy
        real(8) :: temp_x, temp_y, temp_vx, temp_vy
        
        ! k1 coefficients
        call calculate_derivatives(x, y, vx, vy, k1x, k1y, k1vx, k1vy, drag_factor)
        
        ! k2 coefficients  
        temp_x = x + 0.5d0 * dt * k1x
        temp_y = y + 0.5d0 * dt * k1y
        temp_vx = vx + 0.5d0 * dt * k1vx
        temp_vy = vy + 0.5d0 * dt * k1vy
        call calculate_derivatives(temp_x, temp_y, temp_vx, temp_vy, k2x, k2y, k2vx, k2vy, drag_factor)
        
        ! k3 coefficients
        temp_x = x + 0.5d0 * dt * k2x
        temp_y = y + 0.5d0 * dt * k2y  
        temp_vx = vx + 0.5d0 * dt * k2vx
        temp_vy = vy + 0.5d0 * dt * k2vy
        call calculate_derivatives(temp_x, temp_y, temp_vx, temp_vy, k3x, k3y, k3vx, k3vy, drag_factor)
        
        ! k4 coefficients
        temp_x = x + dt * k3x
        temp_y = y + dt * k3y
        temp_vx = vx + dt * k3vx  
        temp_vy = vy + dt * k3vy
        call calculate_derivatives(temp_x, temp_y, temp_vx, temp_vy, k4x, k4y, k4vx, k4vy, drag_factor)
        
        ! Update state variables
        x = x + (dt/6.0d0) * (k1x + 2.0d0*k2x + 2.0d0*k3x + k4x)
        y = y + (dt/6.0d0) * (k1y + 2.0d0*k2y + 2.0d0*k3y + k4y)
        vx = vx + (dt/6.0d0) * (k1vx + 2.0d0*k2vx + 2.0d0*k3vx + k4vx)
        vy = vy + (dt/6.0d0) * (k1vy + 2.0d0*k2vy + 2.0d0*k3vy + k4vy)
    end subroutine runge_kutta_step
    
    ! Calculate derivatives for the differential equation system
    subroutine calculate_derivatives(x, y, vx, vy, dxdt, dydt, dvxdt, dvydt, drag_factor)
        implicit none
        real(8), intent(in) :: x, y, vx, vy, drag_factor
        real(8), intent(out) :: dxdt, dydt, dvxdt, dvydt
        real(8) :: speed, drag_force_factor
        
        ! Position derivatives are simply the velocities
        dxdt = vx
        dydt = vy
        
        ! Calculate drag force components
        speed = sqrt(vx**2 + vy**2)
        
        if (speed > 0.0d0) then
            drag_force_factor = drag_factor * speed
            dvxdt = -drag_force_factor * vx  ! Drag opposes velocity
            dvydt = -g - drag_force_factor * vy  ! Gravity + drag
        else
            dvxdt = 0.0d0
            dvydt = -g
        end if
    end subroutine calculate_derivatives
    
end program projectile_motion
