function varargout = get_phase_index(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(147,varargin{:});
end
