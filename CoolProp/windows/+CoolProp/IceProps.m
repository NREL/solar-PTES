function varargout = IceProps(varargin)
  [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(456,varargin{:});
end
