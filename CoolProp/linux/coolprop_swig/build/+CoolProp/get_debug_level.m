function varargout = get_debug_level(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(318,varargin{:});
end
