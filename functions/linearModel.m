function [Q, S] = linearModel(S, P, K, Dt)
%
% In this function the linear reservoir model is evaluated by using the
% implicite Euler solution.
%
% Input: S  ... storage vector (nx1) with S(1) = initial storage
%        P  ... precipitation input vector (nx1)
%        K  ... model parameter (scalar)
%        Dt ... time step (scalar)
%
% Output: Q ... discharge vector (nx1)
%         S ... storage vector (nx1)
%

% initialize discharge vector (nx1)
Q = zeros(length(P),1);

% Run implicit Euler solution (derived analitically)
for t=1:length(P)
    % Calculate the discharge at the end of the time step using the
    % implicit method (analytical derivation)
    S(t)=(P(t)*Dt+S(t))/(1+K*Dt);
    % Calculate the discharge over the time step (assumed to be the average
    % discharge over the time step)
    Q(t)=K*S(t);
    if t<length(P)
        S(t+1)=S(t);
    end
end

end