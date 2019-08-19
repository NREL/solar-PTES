function varargout = get_config_as_json_string(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(418,varargin{:});
end
