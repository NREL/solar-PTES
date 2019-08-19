function varargout = saturation_ancillary(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(322,varargin{:});
end
