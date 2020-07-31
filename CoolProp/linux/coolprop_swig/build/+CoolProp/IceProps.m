function varargout = IceProps(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(426,varargin{:});
end
