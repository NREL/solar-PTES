function varargout = Props1SI(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(315,varargin{:});
end
