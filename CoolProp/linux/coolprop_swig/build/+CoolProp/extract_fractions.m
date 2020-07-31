function varargout = extract_fractions(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(331,varargin{:});
end
