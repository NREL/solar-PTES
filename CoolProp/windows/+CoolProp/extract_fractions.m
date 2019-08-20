function varargout = extract_fractions(varargin)
  [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(360,varargin{:});
end
