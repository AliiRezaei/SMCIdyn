classdef ThreeDOF
    properties(Access=public)
        m % mass of links
        l % length of links
        lc % length of center of mass 
        I % inertia moment
        g % gravity
        q % joints angular position
        dq % joints angular velocity
    end
    properties(Access=private)
        % position of joints
        JointPos0 = [0, 0, 0]'
        JointPos1
        JointPos2
        % position of end-effector
        EndEffectorPos
        NumLinks = 3; % manipulator degree of freedom (Number of Links)
    end
    methods
        % constructor
        function obj = ThreeDOF(mass, length, length_c, Inertia, grav, pos, vel)
            obj.m = mass;
            obj.l = length;
            obj.lc = length_c;
            obj.I = Inertia;
            obj.g = grav;
            obj.q = pos;
            obj.dq = vel;
        end
        function D = MassMatrix(obj)
            DH = obj.DenavitHartenberg(); % dh for center of masses
            % calculate mass matrix 
            D = 0;
            for i = 1:3
                J_COM = GetJacobian(DH{i}); % get center of mass jacobian
                D = D + (J_COM(1:3,1:end)' * obj.m(i) * J_COM(1:3,1:end) + ...
                    J_COM(4:6,1:end)' * obj.I(i) * J_COM(4:6,1:end));
            end
            
            D = simplify(D);
        end
        function C = CoriolisAcc(obj)
            D = obj.MassMatrix(); % for calculate C matrix we need to mass matrix
           % ready to calculate corilis and centrifugal acceleration matrix
            C = zeros(obj.NumLinks,obj.NumLinks,'sym');
            for k=1:3
                for j=1:3
                    for i=1:3
                        C(k,j) = C(k,j)+1/2*(diff(D(k,j),obj.q(i)) + ...
                            diff(D(k,i),obj.q(j)) - ...
                            diff(D(i,j),obj.q(k)))*obj.dq(i);
                    end
                end
            end
            C = simplify(C);
        end
        function G = GravityVector(obj)
            P1 = obj.m(1) * obj.g * obj.lc(1);
            P2 = obj.m(2) * obj.g * (obj.l(1) + obj.lc(2) * sin(obj.q(2)));
            P3 = obj.m(3) * obj.g * (obj.l(1) + obj.l(2) * sin(obj.q(2)) + obj.lc(3) * sin(obj.q(2) + obj.q(3)));
            P = P1 + P2 + P3; % potantial energy
            % gravity vector :
            G = zeros(3,1,'sym');
            for k = 1:3
                G(k) = diff(P,obj.q(k));
            end
        end
        function PlotHandle = Motion(obj, NowPos, viewVal, varargin)
            % NowPos : Robot Current Position
            % viewVal : the plane that selected to plot motion
            % varargin : plot options
            
            % set plot view
            CurrentAx = gca;
            if strcmp(viewVal,'3D')
                view(3)
            elseif strcmp(viewVal,'YX')
                view(CurrentAx,[0 90])
            elseif strcmp(viewVal,'ZY')
                view(CurrentAx,[90 0])
            elseif strcmp(viewVal,'ZX')
                view(CurrentAx,[0 0])
            else
               error('View argument is Invalid. you can use : 3D, YX, ZY, ZX'); 
            end
            % this method plot Robot Motion in selected view
            obj = EvaluatePos(obj,NowPos);
            
            % joints position in 3D space :
            X = [obj.JointPos0(1),obj.JointPos1(1),obj.JointPos2(1),obj.EndEffectorPos(1)];
            Y = [obj.JointPos0(2),obj.JointPos1(2),obj.JointPos2(2),obj.EndEffectorPos(2)];
            Z = [obj.JointPos0(3),obj.JointPos1(3),obj.JointPos2(3),obj.EndEffectorPos(3)];
            
            % links position in 3D space :
            Link1X = [X(1),X(2)];
            Link1Y = [Y(1),Y(2)];
            Link1Z = [Z(1),Z(2)];
            
            Link2X = [X(2),X(3)];
            Link2Y = [Y(2),Y(3)];
            Link2Z = [Z(2),Z(3)];
            
            Link3X = [X(3),X(4)];
            Link3Y = [Y(3),Y(4)];
            Link3Z = [Z(3),Z(4)];
            
            Line1 = line(Link1X,Link1Y,Link1Z,varargin{:}); % link1
            Line2 = line(Link2X,Link2Y,Link2Z,varargin{:}); % link2
            Line3 = line(Link3X,Link3Y,Link3Z,varargin{:}); % link3
            PlotHandle = {Line1,Line2,Line3};
            hold on
            % plot trajectory that followed by end-effector
            plot3(X(end),Y(end),Z(end),'r.')
            hold off
        end
        % this method plot sliding surface
        function PlotSlidSurface(~, t, s)
            figure
            hold on
            grid minor
            leg = cell(1,size(s,1));
            for k = 1:size(s,1)
                plot(t,s(k,:),'LineWidth',1.5)
                leg{k} = sprintf('S_%d',k);
            end
            legend(leg{:})
            xlabel('Time(s)')
            ylabel('S(x)')
            title('Sliding Surface')
        end
        % thhis method plot control signal
        function PlotControlSignal(~, t, u)
            figure
            for k = 1:size(u,1)
                subplot(3,1,k)
                plot(t,u(k,:),'LineWidth',1.5)
                leg = sprintf('u_%d',k);
                legend(leg)
                grid minor
                ylabel('u')
            end
            xlabel('Time(s)')
            suptitle('Control Signals (N.m)')
        end
        % this method plot error or derivative of error
        function PlotError(~, t, e, dt, mode)
            % if mode = 'error' --> plot error func
            % if mode = 'derror' --> plot derivative of error func
            if strcmp(mode,'error')
                err = {'x','y','z'};
                tit = 'Position Error';
            elseif strcmp(mode,'derror')
                err = {'dx','dy','dz'};
                tit = 'Velocity Error';
            else
                error('invalid input')
            end
            figure
            for k = 1:size(e,1)
                subplot(3,1,k)
                plot(t,e(k,:),'LineWidth',1.5);
                [MaxVal,MaxIdx] = max(abs(e(k,:)));
                MaxErr = sprintf('|Max error| = %7.4f\n at time %5.2f s',MaxVal,(MaxIdx-1) * dt);
                ylabel(err{k})
                legend(MaxErr)
                grid minor
            end
            xlabel('Time(s)')
            suptitle(tit)
        end
        % this method plot robot states
        % in this case : states = joints pos and vel
        function PlotStates(~, t, state)
            JointsPos = state(:,1:3);
            JointsVel = state(:,4:6);
            figure
            for k = 1:size(JointsPos,2)
                subplot(3,1,k)
                plot(t,JointsPos(:,k),'LineWidth',1.5)
                leg = sprintf('q_%d',k);
                legend(leg)
                grid minor
                ylabel('q')
            end
            xlabel('Time(s)')
            suptitle('Joints Position (rad)')
            figure
            for k = 1:size(JointsVel,2)
                subplot(3,1,k)
                plot(t,JointsVel(:,k),'LineWidth',1.5)
                leg = sprintf('dq_%d',k);
                legend(leg)
                grid minor
                ylabel('dq')
            end
            xlabel('Time(s)')
            suptitle('Joints Velocity (rad/s)')
        end
    end % end of public methods
    methods(Access=private)
        function dh = DenavitHartenberg(obj)
            dh1 = [obj.q(1),obj.lc(1),0,0]; % center of mass 1 dh
            dh2 = [obj.q(1),obj.l(1),0,pi/2 % center of mass 2 dh
                obj.q(2),0,obj.lc(2),0];
            dh3 = [obj.q(1),obj.l(1),0,pi/2 % center of mass 3 dh
                obj.q(2),0,obj.l(2),0
                obj.q(3),0,obj.lc(3),0];
            dh = {dh1,dh2,dh3};
        end % end of DenavitHartenberg
        function obj = EvaluatePos(obj, AngPos)
            % homogeneous transform matrices
            T01 = ForwardKinematics([AngPos(1),obj.l(1),0,pi/2]);
            T12 = ForwardKinematics([AngPos(2),0,obj.l(2),0]);
            T23 = ForwardKinematics([AngPos(3),0,obj.l(3),0]);
            T02 = T01 * T12;
            T03 = T02 * T23;
            % update position of joints
            obj.JointPos1 = T01(1:3,end);
            obj.JointPos2 = T02(1:3,end);
            obj.EndEffectorPos = T03(1:3,end);
        end % end of EvaluatePos
    end % end of private methods
end % end of class