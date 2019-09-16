classdef Configuration < SwigRef
  methods
    function self = Configuration(varargin)
      if nargin~=1 || ~ischar(varargin{1}) || ~strcmp(varargin{1},'_swigCreate')
        % How to get working on C side? Commented out, replaed by hack below
        %self.swigInd = CoolPropMATLAB_wrap(439, varargin{:});
        tmp = CoolPropMATLAB_wrap(439, varargin{:}); % FIXME
        self.swigInd = tmp.swigInd;
        tmp.swigInd = uint64(0);
      end
    end
    function delete(self)
      if self.swigInd
        CoolPropMATLAB_wrap(440, self);
        self.swigInd=uint64(0);
      end
    end
    function varargout = get_item(self,varargin)
      [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(441, self, varargin{:});
    end
    function varargout = add_item(self,varargin)
      [varargout{1:nargout}] = CoolPropMATLAB_wrap(442, self, varargin{:});
    end
    function varargout = get_items(self,varargin)
      [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(443, self, varargin{:});
    end
    function varargout = set_defaults(self,varargin)
      [varargout{1:nargout}] = CoolPropMATLAB_wrap(444, self, varargin{:});
    end
  end
  methods(Static)
  end
end
