classdef CriticalState < CoolProp.SimpleState
  methods
    function varargout = stable(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(137, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(138, self, varargin{1});
      end
    end
    function self = CriticalState(varargin)
      self@CoolProp.SimpleState('_swigCreate');
      if nargin~=1 || ~ischar(varargin{1}) || ~strcmp(varargin{1},'_swigCreate')
        % How to get working on C side? Commented out, replaed by hack below
        %self.swigInd = CoolPropMATLAB_wrap(139, varargin{:});
        tmp = CoolPropMATLAB_wrap(139, varargin{:}); % FIXME
        self.swigInd = tmp.swigInd;
        tmp.swigInd = uint64(0);
      end
    end
    function delete(self)
      if self.swigInd
        CoolPropMATLAB_wrap(140, self);
        self.swigInd=uint64(0);
      end
    end
  end
  methods(Static)
  end
end
