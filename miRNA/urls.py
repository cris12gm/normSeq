from django.contrib import admin
from django.urls import include, re_path
from django.conf import settings
from django.conf.urls.static import static
from .views import miRNAResults,queryPlotHTML

urlpatterns = [
    re_path(r'^ajax_plot$', queryPlotHTML, name='ajax_plot'),
    re_path('', miRNAResults.as_view(), name='mirnaresults'),
    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)