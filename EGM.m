classdef EGM
    methods (Static)
        function [M,C,Lhat] = gridpoints(s)
            % create tracker matrices for each variable. use relationship
            % between C and L derived to create labor vectors.
            C=(0:s.n-1)'; 
            M=(0:s.n-1)';
            L=EGM.init_labor(C,M,s);
            Chat = (0:s.n-1)';
            Mhat = (0:s.n-1)';
            Lhat = EGM.init_labor(Chat,Mhat,s);
            
            for l=1:s.numPds
              % calculate utility from each grid point in AlphaVec, using
              % Euler equation with given values of M, C, and L. "Hat"
              % indicates the vectors under the consumer's expectation that
              % they are rationally allocating resources despite the
              % underlying Beta-Delta discounting process.
              [util, labor] = EGM.marg_val_end_per(s, Mhat, Chat, Lhat);
              % set labor to value calculated by interpolation
              LHatVec = labor;
              % calculate the consumption through the inverse utility
              % relationship from the utility calculated above, under an
              % exponetial discounting paradigm
              ChiHatVec=EGM.inv_util(util,labor,s);
              % calculate the consumption under beta-delta discounting
              ChiVec = EGM.inv_util(s.Beta.*util,labor,s);
              % use alpha grid and chi estimates to create mu vectors from
              % assumped behaviors.
              MuVec = s.AlphaVec+ChiVec-s.W*LHatVec;
              % calculate the ratioanlly adjusted market goods vector
              % using the "hat" vectors
              MuHatVec = s.AlphaVec+ChiHatVec-s.W*LHatVec;
              % calculate the labor supply based on the consumption vector
              % from the inverse utility on the rationally discounted value
              LVec = EGM.calc_labor(ChiVec, s);
              
              % append to original matrices
              M=[M, MuVec'];
              C=[C,ChiVec'];
              Chat = [Chat, ChiHatVec'];
              Mhat=[M, MuHatVec'];
              Lhat=[Lhat, LHatVec'];
              L=[L,LVec'];
              
            end;
        end

        function [V_EUP, labor] =marg_val_end_per(s, M, C, L)
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
                % interpolate for labor
                labor = EGM.next_per_labor(m, M, L);
                % integrate into expected utility
                V_EUP = V_EUP+s.Delta.*s.R.*EGM.MU(consumption, labor,s).*prob;
                
            end;
        end

        function c = inv_util( util, labor, s )
            % calculates the consumption function from a variable denoting
            % the utility of a good
            c = (util./(s.Sigma.*(1-labor).^((1-s.Sigma)* ...
                (1-s.Rho)))).^(1/(s.Sigma*(1-s.Rho)-1));

        end

        function util = MU( c, labor, s )
            % calculate the marginal utility derived by a consumer from a
            % certain level of consumption and a certain level of labor
            util = s.Sigma.*c.^(s.Sigma.*(1-s.Rho)-1).* ...
                (1-labor).^((1-s.Sigma)*(1-s.Rho));

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
        
        function util = MU_l( c, labor, s )
            % calculate the marginal utility derived by a consumer from a
            % certain level of consumption and a certain level of labor
            util = (1-s.Sigma).*(-c).^s.Sigma.*(1-labor).^(-s.Sigma) ...
                        .*(c.^s.Sigma .* (1-labor).^(1-s.Sigma)).^(-s.Rho);
        end
        
        function l = calc_labor(c, s)
           % exploits the relationship between labor and consumption to
           % calculate the value of the labor supply from the level of
           % consumption
           l = 1 - c.*(1-s.Sigma)./(s.W*s.Sigma); 
        end
        
        function l = init_labor(c, m, s)
           % exploits the relationship between labor and consumption to
           % calculate the value of the labor supply from the level of
           % consumption
          l = (c - m)./s.W;
          if EGM.MU_l(c,l,s)/s.W ~= EGM.MU(c,l,s)
              l = (c - m)./s.W;
          end
        end
    end
end