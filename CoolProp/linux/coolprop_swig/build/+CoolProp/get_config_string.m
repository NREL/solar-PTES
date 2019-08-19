function varargout = get_config_string(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(417,varargin{:});
end
