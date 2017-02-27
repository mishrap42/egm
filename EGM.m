classdef EGM
    methods (Static)
        function [M,C,L] = gridpoints(s)
            % create tracker matrices for each variable. use relationship
            % between C and L derived to create labor vectors.
            M=(0:s.n-1)';
            [C,L] = EGM.initialize(M,s);
            Chat = C;
            Lhat = L;
            Mhat = M;
                        
            for l=1:s.numPds
              % calculate utility from each grid point in AlphaVec, using
              % Euler equation with given values of M, C, and L. "Hat"
              % indicates the vectors under the consumer's expectation that
              % they are rationally allocating resources despite the
              % underlying Beta-Delta discounting process.
              util = EGM.marg_val_end_per(s, Mhat, Chat);
              % calculate the consumption through the inverse utility
              % relationship from the utility calculated above, under an
              % exponetial discounting paradigm
              ChiHatVec=EGM.inv_util(util,s);
              % compute labor supply consistent with that level of consumption
              LHatVec = EGM.calc_labor(ChiHatVec,s);
              % calculate the consumption under beta-delta discounting
              ChiVec = EGM.inv_util(s.Beta.*util,s);
              % compute labor supply consistent with beta-delta consumption choice
              LVec = EGM.calc_labor(ChiVec, s);
              % use alpha grid and chi estimates to create mu vectors 
              MuVec = s.AlphaVec + ChiVec + 1 - s.W*LVec;
              % calculate the rationally adjusted market goods vector
              % using the "hat" vectors
              MuHatVec = s.AlphaVec + ChiHatVec + 1 - s.W*LHatVec;              
              
              % append to original matrices
              M=[M, MuVec'];
              C=[C,ChiVec'];
              Chat = [Chat, ChiHatVec'];
              Mhat=[Mhat, MuHatVec'];
              Lhat=[Lhat, LHatVec'];
              L=[L,LVec'];
              
            end
        end

        function V_EUP =marg_val_end_per(s, M, C)
        % This function tabulates the marginal utility/value function at
        % the end of the period according to the Bellman equation. It
        % entails the use of two other functions, the interpolation
        % functions next_per_cons and next_per_labor
        
            %expected utility tracker matrix
            V_EUP = zeros(size(s.AlphaVec)); 
            for ii=1:length(s.zval)
                % craft exogenous variables and the probabilities of the
                % shock associated with each
                exogenous=s.zval(ii);
                prob = s.zprob(ii);

                % market goods this period including any exogenous shock
                % factors and the current wage rate
                m = s.R.*s.AlphaVec+exogenous*s.W;
                
                % interpolate for consumption
                consumption = EGM.next_per_cons(m, M, C);
                % integrate into expected utility
                V_EUP = V_EUP+s.Delta.*s.R.*EGM.MU_c(consumption,s).*prob;
                
            end;
        end

        function c = inv_util( util, s )
            % calculates the consumption function from a variable denoting
            % the utility of a good
            c = util.^(-1/s.Rho);

        end

        function util = MU_c( c, s )
            % calculate the marginal utility derived by a consumer from a
            % certain level of consumption and a certain level of labor
            util = c.^(-s.Rho);

        end

        function c=next_per_cons(m, M, C)
            % next_per_cons is an interpolation algorithm to determine
            % next-period consumption function C_t+1()
            mtp1=M(:,end);  % data for the next-period market goods function
            ctp1=C(:,end);  % data for the next-period consumption function

            c = zeros(size(m));

            % extrapolate above maximal m:
            iAbove = m>=mtp1(end);
            slopeAbove = (ctp1(end)-ctp1(end-1))/(mtp1(end)-mtp1(end-1));
            c(iAbove) = ctp1(end) + (m(iAbove)-mtp1(end))*slopeAbove;

            % extrapolate below minimal m:
            iBelow = m<=mtp1(1);
            slopeBelow = 1;
            c(iBelow) = ctp1(1) + (m(iBelow)-mtp1(1))*slopeBelow;

            % interpolate:
            iInterp = ~(iAbove | iBelow);
            c(iInterp) = interp1(mtp1,ctp1,m(iInterp));  

        end

        function l=next_per_labor(m, M, L)
            % next_per_labor is an interpolation algorithm to determine
            % next-period labor function L_t+1()
            mtp1=M(:,end);  % data for the next-period market goods function
            ltp1=L(:,end);  % data for the next-period labor supply function

            l = zeros(size(m));

            % extrapolate above maximal m:
            iAbove = m>=mtp1(end);
            slopeAbove = (ltp1(end)-ltp1(end-1))/(mtp1(end)-mtp1(end-1));
            l(iAbove) = ltp1(end) + (m(iAbove)-mtp1(end))*slopeAbove;

            % extrapolate below minimal m:
            iBelow = m<=mtp1(1);
            slopeBelow = 1;
            l(iBelow) = ltp1(1) + (m(iBelow)-mtp1(1))*slopeBelow;

            % interpolate:
            iInterp = ~(iAbove | iBelow);
            l(iInterp) = interp1(mtp1,ltp1,m(iInterp));  

        end
        
        function util = MU_l(labor, s )
            % calculate the marginal utility derived by a consumer from a
            % certain level of consumption and a certain level of labor
            util = -labor.^(1/s.Epsilon);
        end
        
        function l = calc_labor(c, s)
           % exploits the relationship between labor and consumption to
           % calculate the value of the labor supply from the level of
           % consumption
           l = (s.W./c.^(s.Rho)).^s.Epsilon; 
        end
        
        function [c,l] = initialize(m, s)
            
            % jointly computes optimal c and l which exhaust m
            bc = @(c,l,m) m + l.*s.W - 1 - c;
            foc = @(c,l) EGM.MU_c(c,s) + EGM.MU_l(l,s)./s.W;
            
            c = nan(s.n,1);
            l = nan(s.n,1);
            
            for idx = 1:s.n
                % solve for x = [c l]
                x0 = [m(idx)+0.5*s.W 0.5];
                xSol = fsolve(@(x) [bc(x(1),x(2),m(idx)) foc(x(1),x(2))], x0);
                c(idx) = xSol(1);
                l(idx) = xSol(2);
            end
            
        end
    end
end