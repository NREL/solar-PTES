function varargout = get_config_double(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(416,varargin{:});
end
