classdef GuessesStructure < SwigRef
  methods
    function varargout = T(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(172, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(173, self, varargin{1});
      end
    end
    function varargout = p(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(174, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(175, self, varargin{1});
      end
    end
    function varargout = rhomolar(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(176, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(177, self, varargin{1});
      end
    end
    function varargout = hmolar(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(178, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(179, self, varargin{1});
      end
    end
    function varargout = smolar(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(180, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(181, self, varargin{1});
      end
    end
    function varargout = rhomolar_liq(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(182, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(183, self, varargin{1});
      end
    end
    function varargout = rhomolar_vap(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(184, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(185, self, varargin{1});
      end
    end
    function varargout = x(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(186, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(187, self, varargin{1});
      end
    end
    function varargout = y(self, varargin)
      narginchk(1, 2)
      if nargin==1
        nargoutchk(0, 1)
        varargout{1} = CoolPropMATLAB_wrap(188, self);
      else
        nargoutchk(0, 0)
        CoolPropMATLAB_wrap(189, self, varargin{1});
      end
    end
    function self = GuessesStructure(varargin)
      if nargin~=1 || ~ischar(varargin{1}) || ~strcmp(varargin{1},'_swigCreate')
        % How to get working on C side? Commented out, replaed by hack below
        %self.swigInd = CoolPropMATLAB_wrap(190, varargin{:});
        tmp = CoolPropMATLAB_wrap(190, varargin{:}); % FIXME
        self.swigInd = tmp.swigInd;
        tmp.swigInd = uint64(0);
      end
    end
    function delete(self)
      if self.swigInd
        CoolPropMATLAB_wrap(191, self);
        self.swigInd=uint64(0);
      end
    end
  end
  methods(Static)
  end
end
