function varargout = HAProps(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(424,varargin{:});
end
