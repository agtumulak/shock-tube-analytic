function [ ] = shocktubeanalytic( p_l, T_l, U_l, p_r, T_r, U_r )
%SHOCKTUBEANALYTIC Plots analytical solution of shock tube problem
%   A long 1-D tube is divided into two equal regions by a thin diaphragm.
%   The left side has pressure p_l (Pa), temperature T_l (K), velocity U_l
%   (m)/(s). The right side has pressure p_r, temperature T_r, and 
%   velocity U_r.

%% Temporary (PRESSURE ON LEFT MUST BE GREATER THAN RIGHT OTHERWISE YOU MESS UP)
p_l = 2;
T_l = 1;
U_l = 0;
p_r = 1;
T_r = 1;
U_r = 0;

%% Physical constants

% Gas constant (J)/(mol K)
R = 8.3144598;

% Distance of interest (m)
L = 1;
epsilon = L / 1e10;

% Time of interest (s)
t = 0.050;

%% Gas-dependent constants

% Heat capacity ratio (1)
gamma = 1.4;

% Molecular mass (g)/(mol)
M = 28.96; % Air

% Specific gas constant (J)/(kg K)
R_s = R / M * 1000; % Air

%% Find mach number of shock wave

% Create nondimentionalized variables
den_l_ = Density_( Density( p_l, T_l ) );
den_r_ = Density_( Density( p_r, T_r ) );
T_l_ = Temperature_( T_l );
T_r_ = Temperature_( T_r );
p_l_ = Pressure_( den_l_, T_l_ );
p_r_ = Pressure_( den_r_, T_r_ );
a_l_ = sqrt( gamma * Temperature_( T_r ) );

% Find the root of the compatibility equation to solve for shock mach number
shock_mach_ = fzero( @CompatibilityEquation, 1 );

%% Step 1: Regions 1 and r

% Solve for quantities in region 1
p_1_ = p_r_ * ( ( 2 * gamma * shock_mach_ ^ 2 ) / ( gamma + 1 ) ...
    - ( gamma - 1 ) / ( gamma + 1 ) );
den_1_ = den_r_ / ( 2 / ( ( gamma + 1 ) * ( shock_mach_ ^ 2 ) ) ...
    + ( gamma - 1 ) / ( gamma + 1 ) );
U_1_ = 2 / ( gamma + 1 ) * ( shock_mach_ - 1 / shock_mach_ );

%% Step 2 and 3: Region 2

% Create nondimentionalized variables
den_l_ = Density_( Density( p_l, T_l ) );

% Solve for quantities in region 2
U_2_ = U_1_;
p_2_ = p_1_;
den_2_ = den_l_ * power( p_2_ / p_l_, 1 / gamma );
a_2_ = a_l_ - U_2_ * ( gamma - 1 ) / 2;

%% Step 4: Define compatibility equation to be solved implicitly for shock mach number

    function output = CompatibilityEquation( shock_mach )
        % Define some arguments
        base = p_r_ / p_l_ * ...
            ( ( 2 * gamma * shock_mach ^ 2 ) / ( gamma + 1 ) - ( gamma - 1 ) / ( gamma + 1 ) );
        exponent = ( gamma - 1 ) / ( 2 * gamma );
        % The actual equation
        output = shock_mach - ( 1 / shock_mach ) ...
            - a_l_ * ( gamma + 1 ) / ( gamma - 1) * ( 1 - power( base, exponent ) );
    end

%% Step 5: Calclate positions x_1_, x_2_, x_3_, x_4_

% Create nondimentionalized variables
t_ = Time_( t );

% Solve for x_1_, x_2_, x_3_, x_4_
x_1_ = - a_l_ * t_
x_2_ = ( U_2_ - a_2_ ) * t_
x_3_ = U_2_ * t_
x_4_ = shock_mach_ * t_

%% Return the values of U_e_, a_e_, and p_e_ in the expansion fan

    % Expansion fan speed
    function speed_ = ExpansionFan_U_e_( x_ )
        speed_ = 2 / ( gamma + 1 ) * ( a_l_ + x_ / t_ );
    end

    % Expansion fan sound speed
    function sound_speed_ = ExpansionFan_a_e_( x_ )
        sound_speed_ = a_l_ - ( gamma - 1 ) * ExpansionFan_U_e_( x_ ) / 2;
    end

    % Expansion fan pressure
    function pressure_ = ExpansionFan_p_e_( x_ )
        pressure_ = p_l_ * power( ExpansionFan_a_e_( x_ ) / a_l_, 2 * gamma / ( gamma - 1 ) );
    end

    % Expansion fan density
    function density_ = ExpansionFan_den_e_( x_ )
        density_ = gamma * ExpansionFan_p_e_( x_ ) ./ power( ExpansionFan_a_e_( x_ ), 2);
    end

%% Plot the whole thing

clf
hold on

% Set domain
still_space = 10;
points_l = [-L x_1_];
points_e = linspace( x_1_, x_2_, 1000 );
points_2 = [ x_2_ x_3_ ];
points_1 = [ x_3_ x_3_+epsilon x_4_ ];
points_r = [x_4_ x_4_+epsilon L];

% Get nondimensional velocity U
U_l_region = [0 0];
U_e_region = ExpansionFan_U_e_( points_e );
U_2_region = [U_2_ U_2_];
U_1_region = [U_1_ U_1_ U_1_];
U_r_region = [U_1_ 0 0];

% Plot nondimensional pressure p
p_l_region = [p_l_ p_l_];
p_e_region = ExpansionFan_p_e_( points_e );
p_2_region = [p_2_ p_2_];
p_1_region = [p_1_ p_1_ p_1_];
p_r_region = [p_1_ p_r_ p_r_];

% Plot nondimensional density den
den_l_region = [den_l_ den_l_];
den_e_region = ExpansionFan_den_e_( points_e );
den_2_region = [den_2_ den_2_];
den_1_region = [den_2_ den_1_ den_1_];
den_r_region = [den_1_ den_r_ den_r_];

% Concat the values
all_points = [points_l points_e points_2 points_1 points_r];
all_U = [U_l_region U_e_region U_2_region U_1_region U_r_region];
all_p = [p_l_region p_e_region p_2_region p_1_region p_r_region];
all_den = [den_l_region den_e_region den_2_region den_1_region den_r_region];

% Plot
plot( all_points, [all_U; all_p; all_den] )

title( sprintf(' Analytic solution at nondimensional time t = %s', t_) );
xlabel( 'Nondimensional position' )
ylabel( 'Nondimensional value' )
axis( [-1 1 0 2] )
grid on
legend( 'Velocity', 'Pressure', 'Density' );

%% Convert physical variables to non-dimensional variables.
% Non-dimensional variables are denoted by a trailing underscore. Reference
% state is taken on right side.

    % Position nondimentionalized
    function position_ = Position_( position )
        position_ = position / L;
    end

    % Time nondimentionalized
    function time_ = Time_( time )
        time_ = time / ( L / SoundSpeed( p_r, T_r ) );
    end

    % Density nondimentionalized
    function density_ = Density_( density )
        density_ = density / Density( p_r, T_r) ;
    end
    
    % Speed nondimentionalized
    function speed_ = Speed_( speed )
        speed_ = speed / SoundSpeed( p_r, T_r );
    end

    % Temperature nondimentionalized
    function temperature_ = Temperature_( temperature )
        temperature_ = temperature / ( gamma * T_r );
    end

    % Pressure nondimentionalized
    function pressure_ = Pressure_( density_, temperature_ )
        pressure_ = density_ * temperature_;
    end

    % Energy nondimentionalized
    function energy_ = Energy_( energy )
        energy_ = energy / ( Density( p_r, T_r ) * SoundSpeed( p_r, T_r ) ^ 2 );
    end

    % Enthalpy nondimentionalized
    function enthalpy_ = Enthalpy_( enthalpy )
        enthalpy_ = enthalpy / ( SoundSpeed( p_r, T_r ) ^ 2 );
    end

%% Calculate dependent quantities

    % Density (kg)/(m^3), Pressure (Pa), Temperature (K)
    function density = Density( pressure, temperature )
        density = pressure / (R_s * temperature );
    end

    % Local speed of sound (m)/(s), Pressure (Pa), Temperature (K)
    function a = SoundSpeed( pressure, temperature )
        a = sqrt( gamma * pressure / Density( pressure, temperature ) );
    end

    % Mach number (1), speed (m)/(s), Pressure (Pa), Temperature (K)
    function mach = MachNumber( speed, pressure, temperature )
        mach = speed / SoundSpeed( pressure, temperature );
    end

    % Total enthalpy (m^2)/(s^2), speed (m)/(s), Pressure (Pa), Temperature
    % (K)
    function enthalpy = Enthalpy( speed, pressure, temperature )
        enthalpy = ( SoundSpeed( pressure, temperature ) ^ 2 ) ...
            / ( gamma - 1 ) + ( 1 / 2 ) * speed ^ 2;
    end

end

