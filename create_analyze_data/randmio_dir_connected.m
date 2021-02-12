function [R,eff] = randmio_dir_connected(R, ITER,myIter, attempts)
%RANDMIO_DIR_CONNECTED    Random graph with preserved in/out degree distribution
%
%   R = randmio_dir_connected(W, ITER);
%   [R eff] = randmio_dir_connected(W, ITER);
%
%   This function randomizes a directed network, while preserving the in-
%   and out-degree distributions. In weighted networks, the function
%   preserves the out-strength but not the in-strength distributions. The
%   function also ensures that the randomized network maintains
%   connectedness, the ability for every node to reach every other node in
%   the network. The input network for this function must be connected.
%
%   Input:      W,      directed (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Olaf Sporns, IU

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Mar 2012: Limit number of rewiring attempts, count number of successful
%             rewirings (Olaf Sporns)


n=size(R,1);
[i,j]=find(R);
K=length(i);
ITER=K*ITER;
ITER = myIter;
% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)));
maxAttempts = attempts;
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        rewire=1;
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break           %all four vertices must be different
            end
        end

        %rewiring condition
        if ~(R(a,d) || R(c,b)) 
            [R(a,d) R(c,b)];
            % Vivek: Does either future wiring already exist? 
            % Vivek: If yes, end, do another attempt, and pick two other edges to swap
            % Vivek: If no, check connectedness.
            % If both conditions do not exist, then we rewire.
            %connectedness condition
            [[R(a,c) R(d,b) R(d,c)]; [R(c,a) R(b,d) R(b,a)];];
            [any([R(a,c) R(d,b) R(d,c)]); any([R(c,a) R(b,d) R(b,a)]);];
            [(any([R(a,c) R(d,b) R(d,c)]) && any([R(c,a) R(b,d) R(b,a)]))];
            ~(any([R(a,c) R(d,b) R(d,c)]) && any([R(c,a) R(b,d) R(b,a)]));
            % when we do get into the next if statment we check to make
            % sure that this edge pair is a good candidate to rewire
            % We do not get into this if statement only if 
            % both of the "any" conditions exist
            % we directly rewire when we do not get into the if statement
            % somehow there is a guarentee of connectedness when one of the elements of 
            % ([R(a,c) R(d,b) R(d,c)]) and one of the elements of 
            % [R(c,a) R(b,d) R(b,a)]) are true
            if ~(any([R(a,c) R(d,b) R(d,c)]) && any([R(c,a) R(b,d) R(b,a)]))            
                % conditions to get into this code:
                % 1. If any of the conditions in one any statement set exist,
                % and none of the conditions in the other any statement
                % exist
                % 2. None of the statements in either any statement exist
                P=R([a c],:);
                
                % save the connectivity information in P, for some
                % conditions testing
                % set the connectivity of P to contain the future
                % connectivity
                P(1,b)=0; 
                P(1,d)=1; % this looks like it is setting R(a,b) = 0, and R(a,d) = 1;
                P(2,d)=0; 
                P(2,b)=1; % this looks like it is setting R(c,d) = 0, and R(c,b) = 1;
                
                % set another variable for some more connectivity testing
                % after doing the switch, here self connections are set
                PN=P;
                PN(1,a)=1; % this is essentially setting R(a,a) = 1
                PN(2,c)=1; % this is essentially setting R(c,c) = 1
                
                % P at this point only contains the future edge swapped
                % connections
                
                while 1 % check the network after the edge swap, to see if it meets the conditions for rewiring
                    % P(1,:)~=0, pulls all the connections going from
                    % neuron a ie P(a,:) lets say neuron 2 and 4. 
                    % R(P(1,:)~=0,:) pulls the connections
                    % from neuron 2 and 4 to all other neurons in the
                    % network.
                    % any(...) tells us if there is a connection from either neurons 2 and 4
                    % to all other neurons in the network
                    
                    P(1,:)=any(R(P(1,:)~=0,:),1);
                    P(2,:)=any(R(P(2,:)~=0,:),1);
                    % now P has a new definition that is
                    % whether the edge swapped network neurons connect to other parts
                    % of the network. P(1,:) indicates for the nodes
                    % activated by P(a,:) what next nodes become activated
                    % in any of the nodes activated by P(a,:).
                    % The ~PN means that we want to analyze the set of
                    % links which do not have a connection from from the
                    % edge swaped source node.
                    % P.*(~PN) will tell us whether the source node of the
                    % edge swap activates a node which in turn activates a node 
                    % not activated by the source node. Effectively, the
                    % new P, that is after P=P.*(~PN) will let us know if
                    % we have a "triple of connection" source node of edge
                    % swap, causing a bunch of nodes activating, and out of
                    % the bunch of nodes activating, one or any of those
                    % activating nodes, activates a node not activated by
                    % the original source node. I will call this last node
                    % activated as the Third node.
                    % ~PN is the things which P is not connected to, that 
                    % is NOT everything in P, along with the
                    % self-connection. 
                    % So we now analyze the things that P is not directly 
                    % connected to. 
                    % We look for the places where the successors of P,
                    % the "second layer" are connected to things that the
                    % "first layer", PN, is not connected. In the next
                    % check we are explicitly looking for connections which
                    % "go somewhere else" away from "layer 1" and "layer 2"
                    % New P which is P=P.*(~PN) is the connections from
                    % "layer 2" originating from "layer 1" original P(1,:),
                    % that goes somewhere other than "layer 2" itself, or
                    % back to "layer 1". The condition of "layer 1" is why
                    % the self connection is included in PN.
                    % The new P, now put succinctly, contains the "layer 3"
                    % connections originating from "layer 1". The first row
                    % of the new P is the "layer 3" connection eminating
                    % from a, and the second row are the "layer 3"
                    % connections eminating from c.
                    P=P.*(~PN);
                    PN=PN+P;
                    % An open question, why do the second check in the if
                    % statement??? any(PN...)
                    % The new definition of PN is the connections from
                    % "layer 1" and "layer 2" to the rest of the network.
                    % to PN(1,:) contains the connections from a to 1. itself
                    % 2. "layer 2" and 3. connections which eminate from "layer 2"
                    % Likewise for PN(2,:) and b.
                    if ~all(any(P,2))
                        % We get into this condition if:
                        % 1. A node is activated by both P(a,:),
                        % and P(c,:) then we do not enter this condition.
                        % If there is no Third node activated by both
                        % P(a,:) and P(c,:) then we do not do this
                        % rewiring, because we will have created a "local"
                        % loop. That is a the new P only contains zeros in
                        % P(1,:), P(2,:) or in both.
                        % so a result of any(P,2)=[0;0], any(P,2)=[1;0], 
                        % and any(P,2)=[0;1] would end us up here, a result 
                        % of P=[1;1] would get us to look at the next condition
                        rewire=0;
                        break
                    elseif any(PN(1,[b c])) && any(PN(2,[d a]))
                        % PN gets the new definition
                        % the first statement any(PN(1,[b c])) asks whether
                        % "layer 1" or "layer 2" is connected to b or
                        % c. Both elements of "layer 1" are connected to (b
                        % OR c) AND (d OR a) respectively, then rewire.
                        % PN(1,[b c]) is asking whether node a and/or the 
                        % down stream activated node is connected
                        % to b or c, if a and/or its downstream node is connected to b and/or c, and
                        % node c is connected to a and/or d, then we can
                        % rewire
                        break
                    end
                    % Open question, what happens if we do not meet either
                    % of the criteria in the IF statement?
                    % 
                end
            end %connectedness testing

            if rewire               %reassign edges
                R(a,d)=R(a,b); R(a,b)=0; 
                R(c,b)=R(c,d); R(c,d)=0;
                
                j(e1) = d;          %reassign edge indices
                j(e2) = b;
                eff = eff+1;
                break;
            end %edge reassignment
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations
