function varargout = PropsSImulti(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(317,varargin{:});
end
