function OneDOF
% Simple system
%
%      |----o----|
%          /
%         /
%        /  -^-theta
%  _____/_________)__
    
    close all
    prop.l = 1;
    prop.mass = 1;
    propr.I = 1;
    
    phi = 0;
    theta = -0.1;
    dTheta = 0;
    dPhi = 0;
    tic
    
    dt = 0.01;

    state = [theta, phi, dTheta, dPhi];


    while(true) 
        dt = toc;
        tic
        state = forwardSimulate(state, prop, dt);
        plotSystem(state, prop);
    end
    
end

function state = forwardSimulate(state, prop, tau, dt)
%Simulates the forward dynamics for a state dt in the future
%using rk4
    
    k1 = stateDerivative(state, prop);
    k2 = stateDerivative(state + (dt/2)*k1, prop);
    k3 = stateDerivative(state + (dt/2)*k2, prop);
    k4 = stateDerivative(state + dt*k3, prop);
    
    state = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    theta = state(1);
    if(theta < -pi)
        theta = theta + 2*pi;
    elseif(theta > pi)
        theta = theta - 2*pi;
    end
    state(1) = theta;
end

function state_dot = stateDerivative(state, prop)
    theta = state(1);
    phi = state(2);
    dTheta = state(3);
    dPhi = state(4);
    
    F = calcForce(prop.mass, theta);
    
    
    if(theta < -pi)
        theta = theta + 2*pi;
    elseif(theta > pi)
        theta = theta - 2*pi;
    end

    ddTheta = norm(F)* sign(theta) / (prop.mass * prop.l);

    state_dot = [dTheta, dPhi, ddTheta, 0];
end




function F = calcForce(mass, theta)
%Calculates force on 1D system
    g = -9.8;
    Fg = mass*g;
    FtMag = abs(Fg)*sin(theta);
    F = [FtMag*sin(theta), FtMag*cos(theta)];
end


function plotSystem(state, prop)
    foot = [0,0];
    com = foot + prop.l*[sin(state(1)), cos(state(1))];
    poleOff = .25*[cos(state(2)), sin(state(2))];

    pole = [com - poleOff; com + poleOff];
    
    clf;
    plot([foot(1), com(1)], [foot(2), com(2)]);
    hold on;
    plot(pole(:,1), pole(:,2), 'k');

    axis([-1,1, -.5, 1.5])
    toc
    drawnow
end