function y = fcn(tau_in, tau_out, dt)
%#codegen
persistent input_shaft lay_shaft main_shaft clutch_shaft barrel fork

if isempty(input_shaft)
    input_shaft     = ElementInit(0.01);
    clutch_shaft    = ElementInit(0.015);
    lay_shaft       = ElementInit(0.05);
    main_shaft      = ElementInit(0.07);
    barrel          = ElementInit(0.07);
    fork            = ElementInit(0.01);
end

% Define properties
r_barrel_fork       = 0.01 / deg2rad(45);


% Define model topology and eom's
input_shaft                 = ApplyTorque(input_shaft, tau_in);

[input_shaft, clutch_shaft] = Clutch(input_shaft, clutch_shaft, damping, tau_max);

[clutch_shaft, lay_shaft]	= SpringDamper(clutch_shaft, lay_shaft, 1000, 5, 2);

% Define barrel and dog clutch topology
[barrel, fork]              = SpringDamper(barrel, fork, 1000, 5, r_barrel_fork);

[lay_shaft, main_shaft, fork, a_dog_rel, state] = DogClutch(lay_shaft, main_shaft, 35000, 50, ratio_ab, num_dogs, a_backlash, z_gap, z_fork_offset, fork);

main_shaft                  = ApplyTorque(main_shaft, tau_out);


% State update
[ input_shaft, tau_input_shaft ]   	= ElementUpdate(input_shaft, dt);
[ clutch_shaft, tau_clutch_shaft ]	= ElementUpdate(clutch_shaft, dt);
[ lay_shaft, tau_lay_shaft ]      	= ElementUpdate(lay_shaft, dt);
[ main_shaft, tau_main_shaft ]      = ElementUpdate(main_shaft, dt);

[ barrel, tau_barrel ]              = ElementUpdate(barrel, dt);
[ fork, tau_fork ]                  = ElementUpdate(main_shaft, dt);


% Output values
y   = [ input_shaft.ag
        lay_shaft.ag
        main_shaft.ag
        clutch_shaft.ag
        input_shaft.vel
        lay_shaft.vel
        main_shaft.vel
        clutch_shaft.vel
        ];

end


function [elem] = ElementInit(inertia)
elem.j      = inertia;
elem.ag     = 0;
elem.vel    = 0;
elem.tau    = 0;
end


function [elem] = ApplyTorque(elem, tau)
elem.tau    = elem.tau + tau;
end


function [elem] = EndStopHard(elem, min_pos, max_pos)
if elem.ag < min_pos
    elem.ag     = min_pos;
    elem.vel    = 0;
%     elem.tau    = 0;
elseif elem.ag > max_pos
    elem.ag     = max_pos;
    elem.vel    = 0;
%     elem.tau    = 0;
end
end


function [elem] = GapStopHard(elem, gap, offset)
if abs(elem.ag - offset) < gap       
	elem.ag     = sign(elem.ag - offset) * gap + offset;
	elem.vel    = 0;
end
end


function [elem, tq_elm] = ElementUpdate(elem, dt)
tq_elm      = elem.tau;                         % torque ouput for logging purposes
dxdt        = [ elem.vel; elem.tau / elem.j ];  % state derivative
elem.ag     = elem.ag   + dxdt(1) * dt;         % state integration
elem.vel    = elem.vel  + dxdt(2) * dt;
elem.tau    = 0;                                % reset internal torque
end


function [elem_a, elem_b] = SpringDamper(elem_a, elem_b, stiffness, damping, ratio_ab)
% Compute torque transfer
tau_int     = (elem_a.ag - elem_b.ag / ratio_ab) * stiffness + (elem_a.vel - elem_b.vel / ratio_ab) * damping;
% Apply the reaction torque to both elements
elem_a      = ApplyTorque(elem_a, -tau_int);
elem_b      = ApplyTorque(elem_b, +tau_int * ratio_ab);
end


function [elem_a, elem_b] = Clutch(elem_a, elem_b, damping, tau_max)
% Compute torque transfer
tau_int     = (elem_a.vel - elem_b.vel) * damping;
tau_int     = min(max(tau_int, -tau_max), tau_max);
% Apply the reaction torque to both elements
elem_a      = ApplyTorque(elem_a, -tau_int);
elem_b      = ApplyTorque(elem_b, +tau_int);
end


function [elem_a, elem_b, fork, a_dog_rel, state] = DogClutch(elem_a, elem_b, stiffness, damping, ratio_ab, num_dogs, a_backlash, fork, z_gap, z_fork_offset)
persistent dog_state n_dog_offset
if isempty(dog_state)
    dog_state       = 0;
    n_dog_offset    = 0;
end

% Calculate dog clutch window
a_dog       = elem_a.ag  - elem_b.ag  / ratio_ab;
v_dog       = elem_a.vel - elem_b.vel / ratio_ab;

% Determine whether the dogs are in the window
i_dog_window  = abs(mod(a_dog, 2 * pi / num_dogs) - pi / num_dogs) <= a_backlash;

% Dog engagement logic
if dog_state == 0 && i_dog_window && abs(fork.ag - z_fork_offset) < z_gap
    % Compute number of relative rotations between the dogs
    n_dog_offset    = floor(a_dog / (2 * pi / num_dogs));
    dog_state       = 1;
end

% Handle dog 2 dog contact by blocking the fork displacement
if dog_state == 0 && i_dog_window == 0
%     [fork]  = EndStopHard(fork, - z_fork_offset - z_gap, - z_fork_offset + z_gap);
    [fork]  = GapStopHard(fork, z_gap, z_fork_offset);
%     if abs(fork.ag - z_fork_offset) < z_gap       
%         fork.ag     = sign(fork.ag - z_fork_offset) * z_gap + z_fork_offset;
%         fork.vel    = 0;
%     end
end


% Compute relative dog window angle (subtract relative rotations)
a_dog_rel	= a_dog - (2 * pi / num_dogs) * n_dog_offset - pi / num_dogs;

% Compute the dog transfer torque
tq_damp 	= 0;
tq_stif     = 0;
if dog_state == 1
    % Compute backlash and elastic dog contact
    tq_stif	= (max(a_dog_rel - a_backlash, 0) + min(a_dog_rel + a_backlash, 0)) * stiffness;
    
    % Conditional damping depending on backlash
    if abs(a_dog_rel) > a_backlash
        tq_damp	=  v_dog * damping * 1;
    else    % Maybe add some damping in the backlash?
        tq_damp	=  v_dog * damping * 0;
    end
    
    % Conditional damping (damp only in direction of impact)
%     if sign(tq_stif) ~= sign(tq_damp)
%         tq_damp	= 0;
%     end
end

% Compute reaction torque at flange a (with possible ratio from a --> b)
tau_int     = (tq_stif + tq_damp);

% Handle dog disengagement
if abs(fork.ag - z_fork_offset) >= z_gap
    dog_state	= 0;
end

% Apply the reaction torque to both elements
elem_a      = ApplyTorque(elem_a, -tau_int);
elem_b      = ApplyTorque(elem_b, +tau_int * ratio_ab);

% Output dog clutch state
state       = dog_state;
end
