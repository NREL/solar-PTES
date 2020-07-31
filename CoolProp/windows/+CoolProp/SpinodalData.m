classdef SpinodalData < SwigRef
  methods
    function varargout = tau(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(164, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(165, self, varargin{1});
      end
    end
    function varargout = delta(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(166, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(167, self, varargin{1});
      end
    end
    function varargout = M1(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(168, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(169, self, varargin{1});
      end
    end
    function self = SpinodalData(varargin)
      if nargin~=1 || ~ischar(varargin{1}) || ~strcmp(varargin{1},'_swigCreate')
        % How to get working on C side? Commented out, replaed by hack below
        %self.swigInd = CoolPropMATLAB_wrap(170, varargin{:});
        tmp = CoolPropMATLAB_wrap(170, varargin{:}); % FIXME
        self.swigInd = tmp.swigInd;
        tmp.swigInd = uint64(0);
      end
    end
    function delete(self)
      if self.swigInd
        CoolPropMATLAB_wrap(171, self);
        self.swigInd=uint64(0);
      end
    end
  end
  methods(Static)
  end
end
